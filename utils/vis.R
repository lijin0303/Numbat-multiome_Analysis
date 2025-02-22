##### numbat output vis #####
plot_1cloneBulk = function(bulk,chroms,ignore_loh=T){
  
    min_LLR = 5; min_depth = 8; exp_limit = 2
    dot_size = 0.8; dot_alpha = 0.5;text_size = 10
    raster = FALSE
    segsize = 1.5
    if(nrow(bulk)>10^6){frac = 0.8}else{frac=1}
    if (!all(c('state_post', 'cnv_state_post') %in% colnames(bulk))) {
        bulk = bulk %>%
            mutate(state_post = state,
                   cnv_state_post = cnv_state)
    }
    
    # filter events by LLR
    bulk = bulk %>% 
      mutate(
        LLR = ifelse(is.na(LLR), 0, LLR),
        cnv_state_post = ifelse(LLR < min_LLR, 'neu', cnv_state_post),
        state_post = ifelse(LLR < min_LLR, 'neu', state_post))

    # mark clonal LOH
    if ('loh' %in% colnames(bulk)) {
        bulk = bulk %>% mutate(state_post = ifelse(loh, 'del', state_post))
    }
    marker = 'POS'
    marker_label = 'Genomic position'
    import::from(stringr,str_detect)
    if(ignore_loh){
      loh_flag <- (bulk$state_post=="loh")|(bulk$cnv_state_post=="loh")
      bulk$cnv_state_post <- ifelse(loh_flag,"amp",bulk$cnv_state_post)
      bulk$state_post <- ifelse(loh_flag,"amp",bulk$state_post)
    }
    # fix retest states
    bulk = bulk %>%
        mutate(
            theta_level = ifelse(str_detect(state_post, '_2'), 2, 1),
            state_post = ifelse(
                cnv_state_post %in% c('amp', 'loh', 'del'),
                ifelse(p_up > 0.5, paste0(cnv_state_post, '_', theta_level, '_', 'up'),
                       paste0(cnv_state_post, '_', theta_level, '_', 'down')),
                state_post
        ))%>%
      mutate(logFC = logFC - mu)%>% 
      filter(CHROM %in% chroms) %>% 
      sample_frac(frac)
    
    import::from(data.table,as.data.table)
    D = bulk %>%
        mutate(logFC = ifelse(logFC > exp_limit | logFC < -exp_limit, NA, logFC)) %>%
        mutate(pBAF = ifelse(DP >= min_depth, pBAF, NA)) %>%
        mutate(pHF = pBAF) %>%
        as.data.table %>%
        data.table::melt(measure.vars = c('logFC', 'pHF'))
    
    p = ggplot(D,
               aes(x = get(marker), y = value, color = state_post),
               na.rm=TRUE)

    gaps = numbat::gaps_hg38 %>% filter(end - start > 1e+06)
    acen = numbat::acen_hg38

    segs_exclude = rbind(gaps, acen) %>%
      filter(CHROM %in% bulk$CHROM) %>% 
        mutate(CHROM = factor(as.integer(CHROM))) %>%
        dplyr::rename(seg_start = start, seg_end = end) 
    
    if (nrow(segs_exclude) > 0) {
            p = p + 
              geom_rect(inherit.aes = FALSE, data = segs_exclude,
                        aes(xmin = seg_start, xmax = seg_end, 
                            ymin = -Inf, ymax = Inf),
                        fill = "gray95")}

    legend_breaks = c("neu", "loh_up", "loh_down",
                      "del_up", "del_down", 
                      "amp_up", "amp_down",
                      "bamp", "bdel")
    p = p + 
      geom_point(
            aes(shape = str_detect(state_post, '_2'), alpha = str_detect(state_post, '_2')),
            size = dot_size,
            na.rm = TRUE,
            show.legend = TRUE) +
        geom_hline(
            data = data.frame(y = c(0,1), variable = 'pHF'),
            aes(yintercept = y),
            size = 0, alpha = 0) +
        suppressWarnings(scale_alpha_discrete(range = c(dot_alpha, 1))) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
        theme_classic() +
        theme(
            panel.spacing.x = unit(0, 'mm'),
            panel.spacing.y = unit(3, 'mm'),
            panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
            strip.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.title = element_text(size = text_size),
            strip.text = element_text(size = text_size),
            axis.title = element_text(size = text_size),
            legend.text = element_text(size = text_size),
            plot.margin = margin(t = 1, r = 0, b = 1, l = 0, 'cm')
        ) +
        facet_grid(variable ~ CHROM, scales = 'free', space = 'free_x') +
        scale_x_continuous(expand = expansion(add = 5)) +
        scale_color_manual(
            values = numbat:::cnv_colors,
            limits = names(numbat:::cnv_colors),
            breaks = legend_breaks,
            labels = numbat:::cnv_labels[legend_breaks],
            na.translate = FALSE
        ) +
        guides(
            color = guide_legend(title = "CNV state", override.aes = aes(size = 3), ncol = 1),
            fill = 'none', alpha = 'none', shape = 'none') +
        xlab(marker) +
        ylab('')+ 
      geom_hline(data = data.frame(y = c(-exp_limit, exp_limit),
                                   variable = 'logFC'),
             aes(yintercept = y),
            size = 0, alpha = 0)


    segs = bulk %>%
        group_by(CHROM, seg, seg_start, seg_end) %>%
        summarise(phi_mle = mean(phi_mle)) %>% 
        mutate(variable = 'logFC') %>%
        filter(log2(phi_mle) < exp_limit)
    start = 'seg_start'
    end = 'seg_end'

    p = p + 
          geom_segment(
            inherit.aes = FALSE,
            data = segs,
            aes(x = get(start), xend = get(end), 
                y = log2(phi_mle), yend = log2(phi_mle)),
            color = 'darkred',
            size = segsize) +
          geom_hline(data = data.frame(variable = 'logFC'), 
                     aes(yintercept = 0),
                     color = 'gray30', linetype = 'dashed')
        
      
        
    segs2 = D%>%
        filter(!is.na(state_post)) %>% 
        filter(variable == 'pHF') %>%
        mutate(state_post2 = gsub("_[12]","",state_post)) %>% 
        group_by(CHROM, seg, seg_start, seg_end, state_post2,variable) %>%
        summarise(meanphf = mean(value,na.rm=T)) %>% 
        filter(!is.nan(meanphf))
    p = p + 
      geom_segment(
        inherit.aes = FALSE,
        data = segs2,
        aes(x = get(start), xend = get(end), 
            y = meanphf, yend = meanphf),
        color = "black",
        size = segsize) +
      geom_hline(data = data.frame(variable = 'logFC'), 
                 aes(yintercept = 0),
                 color = 'gray30', linetype = 'dashed')
        
    p = p + xlab(marker_label)

    if (raster) {
        p = ggrastr::rasterize(p, layers = 'Point', dpi = 300)
    }

    return(p)
}

plot_cloneBulks = function(bulks,chroms,normal1=T,title = TRUE,title_size = 8,...) {

    if (!'sample' %in% colnames(bulks)) {
        bulks$sample = 1
    }

    plot_list = bulks %>%
        split(.$sample) %>%
        lapply(
            function(bulk) {
                

                sample = unique(bulk$sample)
                n_cells = sum(unique(bulk$n_cells))
                if(sample==1 & normal1){bulk$LLR <- NA}
                p = plot_1cloneBulk(bulk,chroms,...) +
                    theme(
                        title = element_text(size = title_size),
                        axis.text.x = element_blank(),
                        axis.title = element_blank(),
                        plot.margin = margin(t = 0, r = 0, b = 0.25, l = 0, 'cm')
                    )

                if (title) {
                    if (is.null(n_cells)) {
                        title_text = sample
                    } else {
                        title_text = glue('{sample} (n={n_cells})')
                    }
                    p = p + ggtitle(title_text)
                }
                return(p)
            }
        )
    return(plot_list)
}

relabel_graph_clone <- function(G,clone,...){
  import::from(igraph,V,E)
  G <- igraph::induced_subgraph(G, vids = names(clone))
  V(G)$label <- dplyr::recode(V(G)$label,...)
  # V(G)$label <-  dplyr::recode(V(G)$label,"17a,3b_amp"="17a,3b,8a")
  E(G)$to_label <- dplyr::recode(E(G)$to_label,...)
  G = numbat:::label_genotype(G)
  clone_GT <- igraph::as_data_frame(G, "vertices") %>% {setNames(.$GT, .$clone)}
  for(i in names(clone_GT)){
    clone[[i]]$members <- setNames(clone_GT[i],NULL)
  }
  return(list("graph"=G,"clone"=clone))
}

plot_phylogeny = function(G, clone = NULL,add_dist=T,xyProp = 4,sizep = 0.1,...) {
    edge_label_size = 4; node_label_size = 6; node_size = 10; arrow_size = 2
    show_clone_size = TRUE; show_distance = TRUE; legend = TRUE
    edge_label = TRUE; node_label = TRUE;
    horizontal = TRUE; pal = NULL
    import::from(igraph,V,"V<-","E<-")
    checked_output <- relabel_graph_clone(G,clone,...)
    G <- checked_output$graph
    clone_post <- checked_output$clone
    cloneD <- map(clone_post,\(x) table(grepl("GEX",x$cells))) %>% 
      bind_rows() %>% 
      set_colnames(c("ATAC","RNA")) %>% 
      mutate(group=1:n())
    # G = label_genotype(G,clone_post)
    if (!is.null(clone_post)) {
        clone_sizes = map(clone_post,\(x) x$size)
        igraph::V(G)$size = unname(clone_sizes[V(G)$clone])
    }else{
       igraph::V(G)$size = 1
        show_clone_size = FALSE
    }

    if (is.null(pal)) {
        pal = c("#E41A1C", "#377EB8", "#4DAF4A", 
                "#984EA3", "#FF7F00", "#FFFF33",
                "#A65628", "#F781BF")
        getPalette = colorRampPalette(pal)
        pal = c('gray', getPalette(length(V(G))))
    }
    import::from(tidygraph,as_tbl_graph,activate)
    G_df = G %>% as_tbl_graph() %>% mutate(clone = factor(clone))
    superclone <- 'superclone' %in% colnames(
      as.data.frame(
        activate(G_df, 'nodes')
        ))
    if (!superclone) {
        G_df = G_df %>% mutate(superclone = clone)
    }
    import::from(stringr,str_split,str_trunc)
    # add edge length
    if (show_distance) {
        if ((!'length' %in% colnames(as.data.frame(activate(G_df, 'edges'))))) {
            G_df = G_df %>% activate(edges) %>%
                mutate(n_mut = unlist(purrr::map(str_split(to_label, ','), length))) %>%
                mutate(length = n_mut)
        }

    } else {
        G_df = G_df %>% activate(edges) %>% mutate(length = 1)
    }

    if (!edge_label) {
        G_df = G_df %>% activate(edges) %>% mutate(to_label = '')
    }
    
    require(ggraph)
    p = G_df %>%
        ggraph(
            layout = 'dendrogram',
            length = length
        ) +
        geom_edge_elbow(
            aes(label = str_trunc(to_label, 20, side = 'center')),
            vjust = -1,
            hjust = 0,
            arrow = arrow(length = unit(arrow_size, "mm")),
            end_cap = circle(4, 'mm'),
            start_cap = circle(4, 'mm'),
            label_size = edge_label_size
        ) +
        theme_void() +
        scale_x_continuous(expand = expansion(0.2)) +
        scale_color_manual(values = pal, limits = force) +
        guides(color = 'none')
    if (horizontal) {
      p = p + coord_flip() + scale_y_reverse(expand = expansion(0.2))
    } else {
      p = p + scale_y_continuous(expand = expansion(0.2))
    }
    if(add_dist){
      make_pie <- function(x, y, size, groups, n, rownum) {
        angles <- c(0, 2*pi * cumsum(n)/sum(n))
        do.call("rbind", Map(function(a1, a2, g) {
          xvals <- c(0, sin(seq(a1, a2, len = 500)) * size, 0) + x
          yvals <- c(0, cos(seq(a1, a2, len = 500)) * size*xyProp, 0) + y
          data.frame(x = xvals, y = yvals, group = g, rownum = rownum)
        }, head(angles, -1), tail(angles, -1), groups))
      }
      pieD <- p$data%>%
        mutate(xcor = x,ycor=y,group=as.integer(clone),size=unlist(size)) %>% 
        select(xcor,ycor,group,size) %>% 
        inner_join(cloneD,by="group") %>% 
        mutate(size = size*sizep/max(size)) %>% 
        mutate(r = row_number()) %>%
        rowwise() %>%
        group_map(~ with(.x, make_pie(xcor, ycor, size, c("ATAC", "RNA"),
                                      c(ATAC,RNA), r))) %>%
        bind_rows()
      p_pie <- p + geom_polygon(data = pieD,
                                 aes(x, y, fill = group, 
                                     group = interaction(group, rownum)))+
        scale_fill_manual(values = c("RNA"="#CF5F23ff",
                                     "ATAC" = "#279391ff"),
                          name="")
    }
    
    if (!legend) {
        p = p + guides(size = 'none')
    }

    size_max = max(unlist(G_df %>% activate(nodes) %>% pull(size)))

    if (show_clone_size) {
        p = p + geom_node_point(
                aes(color = as.factor(superclone), size = unlist(size))
            ) +
            scale_size_continuous(
                range = c(0, node_size),
                limits = c(0, size_max),
                name = 'Cells'
            )
    } else {
        p = p + geom_node_point(
                aes(color = as.factor(superclone)),
                size = node_size
            )
    }

    if (node_label) {
        if (show_clone_size) {
            p = p + geom_node_text(aes(label = clone, size = unlist(size)/2),
                                   show.legend = FALSE)
        } else {
            p = p + geom_node_text(aes(label = clone), size = node_label_size)
        }
    }
    
    if(add_dist){
      return(list(p,p_pie))
    }else{
      return(p)
    }
}

plot_segs = function(segs,t=1,c=2) {
  chrom_CNV <- segs %>% 
    distinct(CHROM,cnv_states) %>% 
    filter(cnv_states!="neu") %>% 
    pull(CHROM) %>% unique()
  chrom_labeller <- function(chr){
        chr[chr %in% c(21, 22)] = ''
        return(chr)
    }

  p <- ggplot(segs) +
    geom_rect(
        aes(xmin = seg_start, xmax = seg_end, ymin = -0.5, ymax = 0.5, fill = cnv_states)
    ) +
    theme_void() +
    facet_grid(~CHROM, space = 'free_x', scales = 'free',
               labeller = labeller(CHROM = chrom_labeller)) +
    theme(
        panel.spacing = unit(1, 'mm'),
        strip.background = element_blank(),
        strip.text = element_text(angle = 0,size = c),
        plot.margin = margin(0, 0, 0, 0),
        legend.position = 'top'
    ) +
    scale_fill_manual(
      values = numbat:::cnv_colors,
      labels = numbat:::cnv_labels,
      name = 'CN states'
    ) +
    ggrepel::geom_text_repel(
      data = segs %>% filter(CHROM %in% chrom_CNV) %>% filter(cnv_states!="neu"),
        aes(x = (seg_start+seg_end)/2, y = -0.6,
            label = stringr::str_remove(seg_cons, '\\d+'),size=t),
        min.segment.length = 0,
        vjust = 1,
        hjust = 0,
        direction = 'x',
        segment.curvature = -0.1,
        segment.ncp = 3,
        segment.angle = 20,
        segment.inflect = TRUE,
        max.overlaps = 5
    ) +
    scale_y_continuous(
        expand = expansion(add = c(0.5, 0))
    ) +
    scale_x_continuous(
        expand = expansion(mult = 0.05)
    ) +
    guides(fill = 'none')
  return(p)
}

plot_phylogeny_single = function(G, clone_post = NULL,xyProp = 4,sizep = 0.1) {
  edge_label_size = 4; node_label_size = 6; node_size = 10; arrow_size = 2
  show_clone_size = TRUE; show_distance = TRUE; legend = TRUE
  edge_label = TRUE; node_label = TRUE;
  horizontal = TRUE; pal = NULL
  add_dist = F
  import::from(igraph,V)
  cloneD <- map(clone_post,\(x) table(grepl("GEX",x$cells))) %>% 
    bind_rows() %>% 
    set_colnames(c("ATAC","RNA")) %>% 
    mutate(group=1:n())
  # G = label_genotype(G,clone_post)
  if (!is.null(clone_post)) {
    clone_sizes = map(clone_post,\(x) x$size)
    igraph::V(G)$size = unname(clone_sizes[V(G)$clone])
  }else{
    igraph::V(G)$size = 1
    show_clone_size = FALSE
  }
  
  if (is.null(pal)) {
    pal = c("#E41A1C", "#377EB8", "#4DAF4A", 
            "#984EA3", "#FF7F00", "#FFFF33",
            "#A65628", "#F781BF")
    getPalette = colorRampPalette(pal)
    pal = c('gray', getPalette(length(V(G))))
  }
  import::from(tidygraph,as_tbl_graph,activate)
  G_df = G %>% as_tbl_graph() %>% mutate(clone = factor(clone))
  superclone <- 'superclone' %in% colnames(
    as.data.frame(
      activate(G_df, 'nodes')
    ))
  if (!superclone) {
    G_df = G_df %>% mutate(superclone = clone)
  }
  import::from(stringr,str_split,str_trunc)
  # add edge length
  if (show_distance) {
    if ((!'length' %in% colnames(as.data.frame(activate(G_df, 'edges'))))) {
      G_df = G_df %>% activate(edges) %>%
        mutate(n_mut = unlist(purrr::map(str_split(to_label, ','), length))) %>%
        mutate(length = n_mut)
    }
    
  } else {
    G_df = G_df %>% activate(edges) %>% mutate(length = 1)
  }
  
  if (!edge_label) {
    G_df = G_df %>% activate(edges) %>% mutate(to_label = '')
  }
  
  require(ggraph)
  p = G_df %>%
    ggraph(
      layout = 'dendrogram',
      length = length
    ) +
    geom_edge_elbow(
      aes(label = str_trunc(to_label, 20, side = 'center')),
      vjust = -1,
      hjust = 0,
      arrow = arrow(length = unit(arrow_size, "mm")),
      end_cap = circle(4, 'mm'),
      start_cap = circle(4, 'mm'),
      label_size = edge_label_size
    ) +
    theme_void() +
    scale_x_continuous(expand = expansion(0.2)) +
    scale_color_manual(values = pal, limits = force) +
    guides(color = 'none')
  if (horizontal) {
    p = p + coord_flip() + scale_y_reverse(expand = expansion(0.2))
  } else {
    p = p + scale_y_continuous(expand = expansion(0.2))
  }
  if(add_dist){
    make_pie <- function(x, y, size, groups, n, rownum) {
      angles <- c(0, 2*pi * cumsum(n)/sum(n))
      do.call("rbind", Map(function(a1, a2, g) {
        xvals <- c(0, sin(seq(a1, a2, len = 500)) * size, 0) + x
        yvals <- c(0, cos(seq(a1, a2, len = 500)) * size*xyProp, 0) + y
        data.frame(x = xvals, y = yvals, group = g, rownum = rownum)
      }, head(angles, -1), tail(angles, -1), groups))
    }
    pieD <- p$data%>%
      mutate(xcor = x,ycor=y,group=as.integer(clone),size=unlist(size)) %>% 
      select(xcor,ycor,group,size) %>% 
      inner_join(cloneD,by="group") %>% 
      mutate(size = size*sizep/max(size)) %>% 
      mutate(r = row_number()) %>%
      rowwise() %>%
      group_map(~ with(.x, make_pie(xcor, ycor, size, c("ATAC", "RNA"),
                                    c(ATAC,RNA), r))) %>%
      bind_rows()
    p_pie <- p + geom_polygon(data = pieD,
                              aes(x, y, fill = group, 
                                  group = interaction(group, rownum)))+
      scale_fill_manual(values = c("RNA"="#CF5F23ff",
                                   "ATAC" = "#279391ff"),
                        name="")
  }
  
  if (!legend) {
    p = p + guides(size = 'none')
  }
  
  size_max = max(unlist(G_df %>% activate(nodes) %>% pull(size)))
  
  if (show_clone_size) {
    p = p + geom_node_point(
      aes(color = as.factor(superclone), size = unlist(size))
    ) +
      scale_size_continuous(
        range = c(0, node_size),
        limits = c(0, size_max),
        name = 'Cells'
      )
  } else {
    p = p + geom_node_point(
      aes(color = as.factor(superclone)),
      size = node_size
    )
  }
  
  if (node_label) {
    if (show_clone_size) {
      p = p + geom_node_text(aes(label = clone, size = unlist(size)/2),
                             show.legend = FALSE)
    } else {
      p = p + geom_node_text(aes(label = clone), size = node_label_size)
    }
  }
  
  if(add_dist){
    return(list(p,p_pie))
  }else{
    return(p)
  }
}

##### WGS data vis #####
require(GenomicRanges)
seg_stats <- function(segs,metricD,t){
  segs_metrics <- makeGRangesFromDataFrame(segs) %>% 
    findOverlaps(makeGRangesFromDataFrame(metricD),type = "any") %>% 
    as.data.frame()
  metricD_L <- split(segs_metrics$subjectHits,segs_metrics$queryHits)
  mean_majority <- function(x,t){
    l <- x[x<=t]
    if(length(l)< (length(x)/2)){l <- x[x>t]}
    return(mean(l))}
  segs%<>%mutate(val_avg = map_dbl(metricD_L,\(l) mean_majority(metricD$val[l],t))[rownames(segs)])
  return(segs)
}
plot_wgs <-  function(cnv,segCov,hets) {
  
  set.seed(114)
  chroms = 1:22
  frac = 0.2
  alpha = 0.4
  dot_color = 'gray'
  raster <- F
  colnames(cnv)[1] <- "CHROM"
  cnv %<>% mutate(aberrant = copyNumber > 2.5 | copyNumber < 1.5)
  
  diploid_chroms = cnv %>% group_by(CHROM) %>%
    summarise(diploid = all(copyNumber < 2.4 & copyNumber > 1.7)) %>%
    filter(diploid) %>%
    pull(CHROM)
  
  baseline = segCov %>% filter(CHROM %in% diploid_chroms) %>% pull(logR) %>% median
  
  adj_segCov = segCov %>% 
    mutate(logR = logR - baseline) %>%
    filter(logR < 1.5 & logR > -1.5)
  
  sampled_plotD <- hets %>% 
    select(CHROM,Start,info=baf) %>% 
    mutate(variable = "BAF") %>% 
    rbind(adj_segCov %>% 
            select(CHROM,Start,info=logR) %>% 
            mutate(variable = "logR")) %>% 
    filter(CHROM %in% chroms) %>% 
    mutate(CHROM = factor(CHROM,levels=chroms)) %>% 
    mutate(variable = factor(variable,levels=c("logR","BAF"))) %>% 
    sample_frac(frac)
  
  p = ggplot(sampled_plotD) +
    geom_point(
      aes(x = Start, y = info),
      size = 0.1, alpha = alpha, pch = 16) +
    geom_point(
      data = data.frame(tumorBAF = c(0,1)),
      aes(x = 1, y = tumorBAF),
      size = 0, alpha = 0
    ) +
    facet_grid(variable~CHROM, space = 'free_x', scale = 'free') +
    scale_color_manual(values = c('black', 'red')) +
    theme_classic() +
    theme(
      panel.spacing.x = unit(0, 'mm'),
      panel.spacing.y = unit(3, 'mm'),
      panel.border = element_rect(size = 0.5, color = 'gray', fill = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_blank()
    ) +
    guides(color = 'none') +
    xlab('Genomic position') +
    ylab('')
  
  # p = p + geom_hline(yintercept = c(2/3,1/3,0.5), linetype = 'dashed', color = 'darkred', size = 0.5)
  segs_exclude = numbat::gaps_hg38 %>% filter(end - start > 1e6) %>% 
    rbind(numbat::acen_hg38)%>%
    dplyr::rename(seg_start = start, seg_end = end)
  
  p = p + 
    geom_rect(
      inherit.aes = F,
      data = segs_exclude, 
      aes(xmin = seg_start,
          xmax = seg_end,
          ymin = -Inf,
          ymax = Inf),
      fill = 'gray95')
  
  
  # segment for BAF
  af_segD <- seg_stats(cnv,
                       hets %>% mutate(End=Start) %>% dplyr::rename(val=baf),
                       0.5)
  p = p + 
    geom_segment(
      inherit.aes = F,
      data = af_segD %>% filter(CHROM %in% chroms) %>% 
        select(CHROM,Start,End,seg_BAF=val_avg) %>% 
        mutate(seg_AAF = 1-seg_BAF) %>% 
        gather(ab,seg_AF,-CHROM,-Start,-End) %>% 
        mutate(variable="BAF")%>% 
        mutate(variable = factor(variable,levels=c("logR","BAF"))),
      aes(x = Start, xend = End, y = seg_AF, yend=seg_AF),
      size = 1,
      color = 'blue')
  
  #segment for logR
  cov_segD <- seg_stats(cnv,
                        adj_segCov %>% dplyr::rename(val=logR),
                        0)
  p = p + geom_segment(
    inherit.aes = F,
    data = cov_segD %>% filter(CHROM %in% chroms) %>% 
      select(CHROM,Start,End,val_avg,aberrant) %>%
      mutate(variable="logR")%>% 
      mutate(variable = factor(variable,levels=c("logR","BAF"))),
    aes(x = Start, xend = End, 
        y = val_avg, yend = val_avg, color = aberrant),
    size = 1
  )
  
  if (raster) {
    p = ggrastr::rasterize(p, layers = 'Point', dpi = 300)
  }
  
  return(p)
  
}

##### evaluation plots #####
cloneAssignPlot <- function(confusionMat,clones2zoom=NULL){
  
  xylab <- colnames(confusionMat)[c(2,3)]
  colnames(confusionMat)[c(2,3)] <- c("mode1","mode2")
  if(!is.null(clones2zoom)){
    dat <- confusionMat %>%
      filter(mode1 %in% clones2zoom) %>% 
      filter(mode2 %in% clones2zoom) %>% 
      mutate_at(c("mode1","mode2"),\(m){forcats::fct_drop(factor(m,levels=clones2zoom))}) %>% 
      group_by(mode1,mode2) %>% summarise(n=n())
  }else{
    dat <- confusionMat %>%
      group_by(mode1,mode2) %>% summarise(n=n())
  }
  
  dat_norm <- dat %>% 
    group_by(mode1) %>% 
    mutate(m1p = n*100/sum(n)) %>% 
    group_by(mode2) %>% 
    mutate(m2p = n*100/sum(n)) 
  p1 <-  ggplot(dat_norm, aes(mode1, mode2)) +
    geom_tile(aes(fill = n)) +
    geom_text(aes(label = round(n, 1))) +
    scale_fill_gradient(low = "white", high = "orange")+
    theme_bw()+
    labs(x=xylab[1],y=xylab[2],title = "")
  p2 <- ggplot(dat_norm, aes(mode1, mode2)) +
    geom_tile(aes(fill = m1p)) +
    geom_text(aes(label = round(m1p, 0))) +
    scale_fill_gradient(low = "white", high = "darkblue")+
    theme_bw()+
    labs(x=xylab[1],y=xylab[2],title = "")
  p3 <- ggplot(dat_norm, aes(mode1, mode2)) +
    geom_tile(aes(fill = m2p)) +
    geom_text(aes(label = round(m2p, 0))) +
    scale_fill_gradient(low = "white", high = "darkgreen")+
    theme_bw()+
    labs(x=xylab[1],y=xylab[2],title = "")
  return(ggpubr::ggarrange(p1,p2,p3,nrow=1,legend = "top"))
}
seg_evalPlot <- function(metric_plotD){
  p1 <- ggplot(metric_plotD %>% arrange(CNV,runID,cnv_cutoff),
         aes(x = recall, y = precision,
             color = runID,group=runID)) +
    facet_wrap(~CNV,scales = "free")+
    geom_point(size = 3) +
    geom_path()+
    labs(
      title = "Precision-Recall Curve",
      x = "Recall",
      y = "Precision",
      color = "Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    )+
    guides(color = guide_legend(nrow = 2)) +
    scale_color_brewer(palette = "Dark2")
  
  p2 <- ggplot(metric_plotD,
         aes(x = cnv_cutoff, y = f1,
             color = runID,group=runID)) +
    geom_point(size = 3,alpha=0.7) +
    geom_line()+
    facet_wrap(~CNV,scales = "free")+
    labs(
      title = "F1 with different cut -offs ",
      x = "CNV cutoff",
      y = "",
      color = ""
    ) +
    theme_minimal()+
    scale_color_brewer(palette = "Dark2")
  ggarrange(p1,p2,nrow=2,ncol=1,common.legend = T)
}
uupset <- function(sets,nSet,xlab="",ylab="",setname=1.3,label=" ",sets_bar_col = c("turquoise4")){
  # main_bar_col <- c("violetred4")
  # sets_bar_col <- c("turquoise4")
  # matrix_col <- c("slateblue4")
  # shade_col <- c("wheat4")
  pacman::p_load(UpSetR,grid) 
  g <- as.ggplot(expression(upset(sets,nsets = nSet,
             nintersects = nSet*2,
             number.angles = 0, point.size = 3, line.size = 1.7,
             mainbar.y.label = xlab,
             sets.x.label = ylab,
             mb.ratio = c(0.6, 0.4),
             text.scale = c(1.2, 1.2, 1, 1, setname, 1.5),
             order.by = "freq",
             sets.bar.color = sets_bar_col),
  grid.text( label,x = 0.7,y = 0.02,
             gp = gpar(fontsize = 10,fontface = 3))
  ))
  return(g)
}
ggVBS <- function(data, x, y,fillvar,legendpos=c(0.85, 0.15),guider=2){
  import::from(ggpubr,rremove,font)
  ggplot(data, aes({{x}}, {{y}},fill={{fillvar}}))+
    geom_violin(alpha=0.2, position = position_dodge(width = .75),size=1,color=NA) +
    geom_boxplot(notch = F,outlier.size=-1,color="black",lwd=1,alpha = 0.2,show.legend = F)+
    ggbeeswarm::geom_quasirandom(shape = 21,size=rel(2.5), dodge.width = .75,
                                 color="black",alpha=1,show.legend = F)+
    theme_bw() +
    rremove("legend.title")+
    theme(
      axis.line = element_line(colour = "black",size=1),
      axis.ticks = element_line(size=1,color="black"),
      axis.text = element_text(color="black"),
      axis.ticks.length=unit(0.2,"cm"),
      strip.text = element_text(size = rel(1.2)),
      strip.background = element_blank(),
      legend.position = legendpos)+
    font("xylab",size=rel(1.4))+  
    font("xy",size=rel(1.2))+ 
    font("xy.text", size = rel(1.2)) +  
    font("legend.text",size = rel(1.3))+
    ylim(c(0.8,1))+
    guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"),
                               ncol = guider))
 
}

theme_clear <-function(){
  theme <- theme(
         panel.background = element_rect(fill='transparent'),
         plot.background = element_rect(fill='transparent', color=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         legend.background = element_rect(fill='transparent'),
         legend.box.background = element_rect(fill='transparent')
       )
  return(theme)
}
ggVBS2 <- function(data, x, y,fillvar,shapevar,legendpos=c(0.85, 0.15),guider=2){
  import::from(ggpubr,rremove,font)
  ggplot(data, aes({{x}}, {{y}},fill={{fillvar}},shape={{shapevar}}))+
    geom_violin(alpha=0.2, position = position_dodge(width = .75),
                size=1,color=NA,
                shape=NA) +
    # geom_boxplot(notch = F,outlier.size=-1,color="black",
    #              lwd=1,alpha = 0.2,show.legend = F)+
    ggbeeswarm::geom_quasirandom(size=rel(2.5), dodge.width = .75,
                                 color="black",alpha=1,show.legend = F)+
    theme_bw() +
    rremove("legend.title")+
    theme(
      axis.line = element_line(colour = "black",size=1),
      axis.ticks = element_line(size=1,color="black"),
      axis.text = element_text(color="black"),
      axis.ticks.length=unit(0.2,"cm"),
      strip.text = element_text(size = rel(1.2)),
      strip.background = element_blank(),
      legend.position = legendpos)+
    font("xylab",size=rel(1.4))+  
    font("xy",size=rel(1.2))+ 
    font("xy.text", size = rel(1.2)) +  
    font("legend.text",size = rel(1))+
    scale_shape_manual(labels=mode_label,values=mode_shapes)+
    ylim(c(0.8,1))+
    guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"),
                               ncol = guider),
           shape = guide_legend(override.aes = list(alpha = 1,color="black"),
                               ncol = guider))
  
}

