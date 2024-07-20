

#############################################################
# helper functions 											#
#############################################################


create_directories <- function(directory){
# create results folder within directory with date in folder name and output list with corresponding PATHs
	folder <- c('plots','results')
	paths = c()
	for (i in 1:length(folder)) {
		temp <- file.path(directory, paste(folder[i],Sys.Date(),sep=''))
		if (!file.exists(temp)) {
			dir.create(temp)
		}
		paths$folder[i] = temp
	}
	return(paths)
}

set_defaults <- function(split = 'expe.day', identifier = 'id', grouping = 'genotype') {
	defaults <- list()
	defaults$grouping <- grouping
	defaults$split <- split
	defaults$identifier <- identifier
	# my.colours <- cm.colors(length(table(unique(scoring$data[defaults$grouping]))), alpha = .5)
	# 	names(my.colours) <- as.matrix(unique(scoring$data[defaults$grouping]))
	# defaults$colours <- my.colours
	defaults$theme <-theme(legend.title=element_blank(),legend.key = element_blank(), legend.position='none', strip.text.x = element_text(size=10), strip.background = element_rect(colour="black", fill="white"),axis.text.x = element_text(angle=0), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"),panel.background = element_blank(),plot.title = element_text(hjust = 0), axis.title.x=element_blank())
	log_msg(LOG_FILE,'DEFAULTS: initialize global defaults')
	return(defaults)
}

log_msg <- function(filename, message) {
	print(message)
	sink(filename,append=TRUE)
	print(message)
	sink()
}

check_libraries <- function(req_libraries) {
	packages<-installed.packages()[,1]
	for (pkg in req_libraries) {
		if (!is.element(pkg, packages)) {
			install.packages(pkg,repos="https://cran.uni-muenster.de") # default german CRAN mirror
			log_msg(LOG_FILE,paste('LIBRARIES:','installing package: ',pkg, sep='   '))
		}
	}
	log_msg(LOG_FILE,paste('LIBRARIES:','all required libraries installed', sep='   '))
	#print(paste('STATUS:','all required libraries installed', sep='   '))
}

load_libraries <- function(req_libraries) {
	for (pkg in req_libraries) {
		#suppressMessages(eval(bquote(library(.(pkg)))))
		suppressMessages(require(pkg, character.only=TRUE))
	}
	log_msg(LOG_FILE,paste('LIBRARIES:','all libraries loaded into workspace', sep='   '))
}


####################################################
# DATA IMPORT AND MANIPULATIONS					   #
####################################################

get_data <- function(file = NULL) {
	if (is.null(file) | !is.element(file,dir(paste(ORIG,'/Data',sep='')))) {
		file <- tk_select.list(dir(paste(ORIG,'/Data',sep = '')), preselect = NULL, multiple = FALSE, title = 'Select Scoring CSV file')
	}
	scoring = list()
	scoring$file <- paste(ORIG,'/Data/',file,sep = '')
	scoring$data <- read.table(scoring$file, header = TRUE, sep = ';',stringsAsFactors=FALSE, check.names=FALSE)
	if (any(grepl('unknown',scoring$data[[DEFAULTS$identifier]]))){
		scoring$data <- scoring$data[-which(scoring$data[[DEFAULTS$identifier]]=='unknown'),]
	}
	my.colours <- cm.colors(length(table(unique(scoring$data[DEFAULTS$grouping]))), alpha = .5)
	names(my.colours) <- as.matrix(unique(scoring$data[DEFAULTS$grouping]))
	scoring$colours <- my.colours
	log_msg(LOG_FILE,paste('LOADING DATA: ',as.character(file)))
	return(scoring)
}

# mapping of scored behaviours to defined groups as defined in: './Data/test_behaviour.csv'
behave_mapping <- function(scoring) {
	x <- scoring$data
	behav_table <- suppressWarnings(read.table(paste(ORIG,'/Data/test_behaviour.csv',sep=''),sep = ';', header = TRUE, check.names=FALSE,stringsAsFactors=FALSE))
	for (beh in names(behav_table)) {
		x$behaviour[x$behaviour %in% behav_table[[beh]]] <- beh
	}
	scoring$mapped.data <- x[x$behaviour %in% names(behav_table),]
	scoring$behaviours <- names(behav_table)
	return(scoring)
}
# adjust of identifying color code for animals in scored data
color_mapping <- function(x) {
	batches <- get_batches(x)
	for (b in batches) {
		temp = x[x$batch == b,]
		batch.colors <- unique(temp$color)
		for (c in batch.colors) {
			animal.nr <- temp[[DEFAULTS$identifier]][temp$color == c][1]
			temp$with.who[temp$with.who == c] <- animal.nr
		}
		x[x$batch==b,] <- temp
	}
	return(x)
}
# checking match of identifiers color & id in data
check_encoding <- function(scoring) {
	x <- scoring$mapped.data
	if (class(x$color) != class(x[[DEFAULTS$identifier]])) {
		x <- color_mapping(x)
	}
	scoring$mapped.data <- x
	return(scoring)
}


#############################################################
# helper functions for analysis steps and data handling		#
#############################################################

get_batches <- function(data) {
	batches <- unique(data$batch)
	return(batches)
}

get_genotypes <- function(data) {
	genotype <- unique(data$genotype)
	return(genoytpe)
}

get_gt_batch_matching <- function(scoring) {
	batches <- get_batches(scoring)
	gt_batch <- matrix(,ncol = 2, nrow=length(batches))
	colnames(gt_batch) <- c('batch','genotype')
	for (b in batches) {
		gt_batch[which(batches==b),1] = b
		gt_batch[which(batches==b),2] = scoring$genotype[scoring$batch==b][1]
	}
	gt_batch <- as.data.frame(gt_batch[order(gt_batch[,'genotype']),])
	return(gt_batch)
}

# get_behaviours <- function(data) {
# 	behaviours <- unique(data$behaviour)
# 	return(behaviours)
# }

behav_split <- function(temp, behaviour) {
	temp.behav <- temp[temp$behaviour == behaviour,]
	return(temp.behav)
}

batch_split <- function(data, batch) {
	temp <- data[data$batch == batch,]
	return(temp)
}

data_split <- function(temp, split = NULL) {
	if (is.null(split)){
		split <- 'expe.day'
	}
	split.instances <- unique(temp[split])
	data <- list()
	data$split = list()
	data$instances <- as.vector(split.instances)
	for (s in split.instances[[split]]) {
		data$split[[paste(split,s, sep = '_')]] <- temp[temp[split] == s,]
	}
	return(data)
}

get_animals <- function(scoring) {
	animals <- list()
	for (b in scoring$batches) {
		temp <- unique(scoring$mapped.data[[DEFAULTS$identifier]][scoring$mapped.data$batch == b])
		id.unknown <- grepl('unknown', temp)
		if (any(id.unknown)) {
			animals[[b]] <- temp[-which(id.unknown)]
		} else {
			animals[[b]] <- temp
		}
		
	}
	return(animals)
}

colMax <- function(x) sapply(x, max, na.rm = TRUE)
colMedian <- function(x) sapply(x, median, na.rm = TRUE)


remove_empty_behaviour <- function(data) {
	empty.id <- which(data$behaviour == '')
	if (length(empty.id)) {
		data = data[-empty.id,]
	}
	return(data)	
}


#############################################################
# wrapper helper function  for computation of 				#
# social networks and metrics								#
#############################################################

conflict_matrix <- function(x, animals.batch) {
	#create empty matrix for cummulated conflict result values
	size = length(as.matrix(animals.batch[[1]]))
	animals.active <- unique(x[[DEFAULTS$identifier]])
	con <- matrix(0, nrow=size, ncol = size)
	colnames(con) <- as.matrix(animals.batch[[1]])
	rownames(con) <- as.matrix(animals.batch[[1]])
	# compute conflict results for each pair of animals
	for (i in animals.active) {
		temp <- x[x[[DEFAULTS$identifier]]==i,]#filter(x,id==i)
		count <- table(temp$with.who)
		try(con[as.character(i),names(count)] <- count, silent = TRUE)
	rm(temp)
	}
	# delete self directed interactions
	diag(con) <- 0
	# return conflict matrix
	return(con)
}


create_networks <- function(data, animals) {
	networks <- list()
	for (s in 1:nrow(data$instances)) {
		networks[[paste(names(data$instances),data$instances[s,], sep = '_')]] <- conflict_matrix(data$split[[s]],animals)
	}
	return(networks)
}

get_graph_obj <- function(adjacency) {
	graph <-  graph_from_adjacency_matrix(adjacency,weighted=TRUE)
	return(graph)
}

#--------------------graph parameter------------------#
strength_in <- function(graph) {
	str <- mean(strength(graph, mode = 'in'))
	names(str) <- 'In strength'
	return(str)
}

strength_out <- function(graph) {
	str <- mean(strength(graph, mode = 'out'))
	names(str) <- c('Out strength')
	return(str)
}

degree_in <- function(graph) {
	deg<- mean(degree(graph, mode = 'in'))
	names(deg) <- 'In degree'
	return(deg)
}

degree_out <- function(graph) {
	deg<- mean(degree(graph, mode = 'out'))
	names(deg) <- 'Out degree'
	return(deg)
}

deg_centr <- function(graph) {
	deg <- centr_degree(graph, normalized = TRUE)$centralization
	return(deg)
}

short_path <- function(graph) {
	path <- mean_distance(graph, directed = TRUE, unconnected = TRUE)
	return(path)
}

betweenness_avg <- function(graph){
	#defined as: 
	# ‚sum( g_ivj / g_ij, i!=j,i!=v,j!=v)
	bet <- betweenness(graph)
	bet.avg <- mean(bet)
	return(bet.avg)
}

dens <- function(graph) {
	d <- edge_density(graph, loops = FALSE)
	return(d)
}

recip <- function(graph) {
	#reciprocity defines the proportion of mutual connections
	r <- reciprocity(graph, ignore.loops = TRUE, mode = 'default')
	return(r)
}


compute_parameter <- function(data,parameter) {
	observations <- names(data)
	param <- as.data.frame(matrix(,ncol = length(parameter)+1,nrow = length(observations)))
	colnames(param) <- c(DEFAULTS$split,parameter)
	#rownames(param) <- observations
	for (o in observations) {
		graph <- get_graph_obj(data[[o]])
		idx = which(observations==o)
		param[idx,DEFAULTS$split] <- o
		param[idx,'strength.in'] <- strength_in(graph)
		param[idx,'strength.out'] <- strength_out(graph)
		param[idx,'degree.in'] <- degree_in(graph)
		param[idx,'degree.out'] <- degree_out(graph)
		param[idx,'short.path'] <- short_path(graph)
		param[idx,'betweenness'] <- betweenness_avg(graph)
		param[idx,'deg.centr'] <- deg_centr(graph)
		param[idx,'density'] <- dens(graph)
		param[idx,'reciprocity'] <- recip(graph)
	}
	return(param)
}

get_parameter <- function(data,gt_batch) {
	batches <- names(data)
	parameter <- c('strength.in','strength.out','degree.in','degree.out','short.path','betweenness','deg.centr','density','reciprocity')
	par <- as.data.frame(matrix(,ncol = length(parameter)+2,nrow=0))
	colnames(par) <- c('batch',DEFAULTS$split,parameter)
	for (b in batches) {
		observations <- names(data[[b]])
		if (is.null(observations)) {
			temp <- cbind(as.matrix(b),matrix(NA,nrow=1,ncol=length(parameter)+1))
			colnames(temp) <- c('batch',DEFAULTS$split,parameter)
		} else {
			temp <- cbind(rep(b,length(observations)),compute_parameter(data[[b]], parameter))
			colnames(temp) <- c('batch',DEFAULTS$split,parameter)
		}
		# temp <- cbind(rep(as.integer(substr(b,2,3)),length(observations)),compute_parameter(data[[b]], parameter))
		par <- rbind(par,temp)
	}
	# colnames(par) <- c('batch',DEFAULTS$split,parameter)
	par <- par[get_ordering(par,gt_batch),]
	return(par)
}

replace_null_na <- function(matrix,value) {
	matrix[matrix=='NULL'] <- value
	matrix[is.na(matrix)] <- value
	return(matrix)
}

get_ordering <- function(data, gt_batch) {
	batch.order <- rep(1,nrow(data))
	idx <- match(gt_batch$batch, data$batch)
	id.order <- sort(idx,decreasing=FALSE)
	observations <- as.matrix(unique(data[DEFAULTS$split]))
	for (i in idx) {
		batch.order[id.order[which(idx==i)]] = i
		for (o in seq(1,length(observations)-1)) {
			batch.order[id.order[which(idx==i+o)]] = i+o
		}
	}
	return(batch.order)
}

#--------------------graph plotting------------------#

graph_plot <- function(graph, title = NULL, scheme = NULL) {
	#adjust display options
	V(graph)$color = adjustcolor('black', alpha.f = .25)
	V(graph)$label.color = 'white'
	V(graph)$label.cex = 1.25
	V(graph)$frame.color = NA
	V(graph)$size = 25
	E(graph)$arrow.size = .5#(E(g)$weight)/25
	E(graph)$curved = 0.1
	if (is.null(scheme)) {
			layout = layout_in_circle(graph)
	}
	plot(graph, layout = layout)

}

# graph_layout <- function(graph) {
# 	size = length(V(graph))+1
# 	theta <- seq(0,2*pi, length = size)
# 	layout = cbind(cos(theta[-1]),sin(theta[-1]))
# 	return(layout)
# }

# plot_sna_graphs <- function(data) {
# 	behaviours <- names(data)
# 	batches <- names(data[[1]])
# 	observations <- names(data[[1]][[1]])
# 	for (beh in behaviours) {
# 		par(mfrow=c(2,ceiling(length(observations)/2)))
# 		for (b in batches) {
# 			temp <- data[[beh]][[b]]
# 			for (o in observations) {
# 				g <- get_graph_obj(temp[[o]])
# 				graph_plot(g)
# 			}


# 		}

# 	}
# }


binarize <- function(matrix) {
	matrix <- ifelse(matrix > 0,1,0)
	return(matrix)
}


  



circle_plot <- function(data, title, leg.pos, nr.obs) {
	empty_bar <- 4
	to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
	colnames(to_add) <- colnames(data)
	to_add$group <- rep(levels(data$group), each=empty_bar)
	data <- rbind(data, to_add)
	# data <- data %>% arrange(batch)
	data$id <- seq(1, nrow(data))
	 
	# Get the name and the y position of each label
	label_data <- data
	number_of_bar <- nrow(label_data)
	angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
	label_data$hjust <- ifelse( angle < -90, 1, 0)
	label_data$angle <- ifelse(angle < -90, angle+180, angle)
 	bar.colors <- get_circle_colors(data,gt_batch)
	# Make the plot
	p <- ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
	  geom_bar(stat="identity", alpha=0.75) +
	  # scale_colour_gradient(  low = scoring$colours[1],
  		#high = scoring$colours[2]) +
	  scale_fill_manual(values=bar.colors) +
	  # ADD: create color scheme including genotype information by rep(color(KO),rep) using cm.color()
	  ylim(-max(data$value)*2,max(data$value)) +
	  theme_minimal() +
	  theme(
	    legend.position = leg.pos,
	    axis.text = element_blank(),
	    axis.title = element_blank(),
	    panel.grid = element_blank(),
	    plot.margin = unit(rep(-1,4), "cm") ,
	    plot.title = element_text(size = 15)
	  ) +
	  coord_polar() + 
	  geom_text(data=label_data, aes(x=id, y=value, label=rownames(data), hjust=hjust), color="darkgray", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
	  labs(title=title)
	 return(p)

}

circle_wrapper <- function(data, beh) {
	plots <- list()
	observations <- unique(data[[DEFAULTS$split]])
	parameter <- colnames(data)[3:ncol(data)]
	data <- replace_null_na(data,0)
	legend_position = 'none'
	nr.obs = length(observations)
	# par(mfrow=c(3,3))
	for (p in parameter) {
		temp <- data[,c('batch',p)]
		temp <- as.data.frame(cbind(as.factor(rownames(temp)),temp))
		colnames(temp) <- c('observation','group','value')
		temp$value <- 1+temp$value
		if (p == parameter[length(parameter)]) {
			legend_position = 'right'
		}
		plots[[p]] <- circle_plot(temp, title=as.character(p), leg.pos=legend_position,nr.obs)
	}
	grid.arrange(grobs=plots, nrow=2, top = grid::textGrob(paste('Network parameter for',beh), gp=gpar(fontsize=20,font=3,hjust=.5)))
}

get_circle_colors <- function(data,gt_batch){
 	group.size <- table(gt_batch[[DEFAULTS$grouping]])
	color_long = cm.colors(length(gt_batch$batch)*1.5)
	color1 <- rep(color_long[seq(1,group.size[1])],each=length(unique(data[[DEFAULTS$split]])))
	color2 <- rep(color_long[seq(length(color_long)-group.size[2], length(color_long))],each=length(unique(data[[DEFAULTS$split]])))
	bar.colors <- c(color1,color2)
	return(bar.colors)
 }

#############################################################
# wrapper and helper function  for computation of 			#
# glicko glicko ratings										#
#############################################################


get_glicko <- function(data) {
	data <- data[,c('time_start',DEFAULTS$identifier,'with.who','behaviour')]
	data$time_start <- as.numeric(chron::times(data$time_start))
	data$behaviour <- 1
	glicko.rating <- glicko(data, init = c(0,1), cval = 1,
		sort = TRUE, history = TRUE)
	return(glicko.rating)
}

create_glicko_ratings <- function(data) {
	data$behaviour <- 1
	#data <- data[which(data$with.who %in% unique(data[[DEFAULTS$identifier]])),]
	glicko <- get_glicko(data)
	return(glicko)
}

check_interactions <- function(data) {
	if (all(unique(data$with.who) == '') | nrow(data) == 1) {
		glicko <- paste('no dyadic interaction for given behaviours')
	} else {
		glicko <- create_glicko_ratings(data)
	}
	return(glicko)
}

#--------------------glicko plotting------------------#



replace_rankings <- function(data,rankings) {
	data$variable <- as.character(data$variable)
	rankings$variable <- as.character(rankings$variable)
	for (r in rankings$variable) {
		data$variable[data$variable==r] <- which(rankings$variable==r)
	}
	return(data)
}

adjust_data_format <- function(data) {
	if (is.character(data)) {
		return(data)
	} else {
		history = as.data.frame(t(data$history[,,1]))
		games=as.data.frame(seq(1,ncol(data$history[,,1])))
		data_extract=cbind(history,games)
		colnames(data_extract) = c(colnames(history),'games')
		data_form <- melt(data_extract, id.var='games')
		rankings = as.data.frame(data$ratings$Player)
		colnames(rankings)<-c('variable')
		data_new=left_join(data.frame(name=rankings),data_form,by="variable")
		data_new <- replace_rankings(data_new, rankings)
		return(data_new)
	}
}
glicko_ggplot <- function(temp,plottitle,ylabel=FALSE,first='') {
	if (ylabel!=FALSE) {
		ytext = ylabel
	} else {ytext = ''}
	p <- ggplot(data=temp, aes(x=games, y=value, linetype = variable)) +
		geom_line() +
		#ggtitle(title) +
		DEFAULTS$theme +
		labs(title = paste(first),
       subtitle = paste(plottitle),
       y=ytext)# 
	return(p)
}

ggplot_empty <- function(temp, plottitle, ylabel=FALSE,first='') {
	if (ylabel!=FALSE) {
		ytext = ylabel
	} else {ytext = ''}
	p <- ggplot() +
		geom_text(aes(x=1,y=1),label='no interactions') +
		DEFAULTS$theme +
		labs(title = paste(first),
        subtitle = paste(plottitle),
        y=ytext)
    return(p)
}

glicko_ggplot_wrapper <- function(glicko,scoring) {
	gt_match <- get_gt_batch_matching(scoring$data)
	batches <- gt_match$batch
	plots <- list()
	new_row = c(1,(length(batches)/2)+1)
	gtsigns <- c(rep('+/+',(length(batches)/2)),rep('-/-',(length(batches)/2)))

	for (b in batches){
		temp <-adjust_data_format(glicko[[b]])
		idx = which(batches==b)
		if (idx %in% new_row) {
			ylabel='glicko rating'
			# ylabel=gtsigns[idx]
		} else {ylabel=FALSE}
		if (idx == 1) {
			first ='A:Glicko rating history'
		} else {
			first=''
		}
		if (is.character(temp)) {
			p <- ggplot_empty(temp,as.character(idx),ylabel,first)
		} else {
			p <- glicko_ggplot(temp, as.character(idx), ylabel, first)
		}
		plots[[b]] <- p
	}
	return(plots)
}



#################
# glicko ranking correlation
#################


rank_cor_int <- function(data,percent,m) {
	if (grepl('S',m)) {
		method = 'spearman'
	} else {
		method = 'pearson'
	}
	if (is.character(data)) {
		return(data)
	} else {
		data=data$history
		size = dim(data)
		nr.intervals = 1/percent
		if (nr.intervals > size[2]) {
			nr.intervals=size[2]-1
		}
		CM = matrix(,nrow = size[2]-1,ncol = 1)
		for (t in seq(1,size[2]-1)) {
			if (sum(data[,t,1])!=0) {
				CM[t] = cor(data[,t,1],data[,size[2],1], method = method)
			} else {
				CM[t] = 0
			}
				
		}
		istop = ceiling(length(CM)*percent)
		ismin = floor(length(CM)*percent)
		sizedif = length(CM)-(nr.intervals*ismin)
		id = sizedif * istop
		em = matrix(,ncol = nr.intervals,nrow= istop)
		em[1:id] = CM[1:id]
		if (id < nr.intervals) {
			em[seq(1,ismin),seq(sizedif+1,nr.intervals)] = CM[seq(id+1,length(CM))]
		}	
		return(em)
	}
}


glicko_cor_plot <- function(data,plottitle,ylabel=FALSE,first='') {
	if (ylabel!=FALSE) {
		ytext = ylabel
	} else {ytext = ''}
	if (is.character(data)) {
		p <- ggplot_empty(data,plottitle,ylabel,first)
	} else if (length(unique(as.vector(data))) <= 3) {
		data_new <- melt(data)
		p <- ggplot(data = data_new, aes(x = Var2, y=value, group = Var2)) +
		geom_boxplot(na.rm=TRUE) +
		ggtitle(title) +
		DEFAULTS$theme +
		labs(title = paste(first),
      	 subtitle = paste(plottitle),
      	 y=ytext)
	} else {
		data_new <- melt(data)
		p <- ggplot(data = data_new, aes(x = Var2, y=value, group = Var2)) +
		geom_boxplot(na.rm=TRUE) + 
		geom_line(data = as.data.frame(loess.smooth(x=matrix(data_new$Var2), y=rowMeans(matrix(data_new$value),na.rm = TRUE))),aes(x=x,y=y), inherit.aes=FALSE) +
		ggtitle(title) +
		DEFAULTS$theme +
		labs(title = paste(first),
  	     subtitle = paste(plottitle),
   	    y=ytext)
	}
	return(p)
}

rank_cor_wrapper <- function(glicko,percent,scoring) {
	gt_match <- get_gt_batch_matching(scoring$data)
	batches <- gt_match$batch

	output <- list()
	rank_correlations <- list()
	plots <- list()
	new_row = c(1,(length(batches)/2)+1)
	gtsigns <- c(rep('+/+',(length(batches)/2)),rep('-/-',(length(batches)/2)))

	for (b in batches){
		idx = which(batches==b)
		if (idx %in% new_row) {
			ylabel='rank correlation'
		} else {ylabel=FALSE}
		if (idx == 1) {
			first ='B:Rank correlations'
		} else {
			first=''
		}
		rank_correlations[[b]] <- rank_cor_int(glicko[[b]],percent,'spearman')

		plots[[b]] <- glicko_cor_plot(rank_correlations[[b]],as.character(idx),ylabel,first)
	}
	output[['correalations']] <- rank_correlations
	output[['plots']] <- plots
	return(output)
}



arange_plots <- function(plots, title, yaxistitle, xaxistitle) {
	grid.arrange(grobs=plots, nrow=2, 
		top = grid::textGrob(as.character(title),gp=gpar(fontsize=20,font=3,hjust=0.5)), 
		bottom = grid::textGrob(xaxistitle,gp=gpar(fontsize=10,font=3,hjust=0)),
		left = grid::textGrob(as.character(yaxistitle),gp=gpar(fontsize=10,font=3,hjust=0,rotate=90)))
}

#################
# glicko despotism plotting
#################

plot_despotism <- function(data) {
	groups <- unique(data[,DEFAULTS$grouping])
	nr_int <- 100 # set default value for interpolation bins in density estimation
	power <- matrix(,nrow = length(groups),ncol = nr_int)
	rownames(power) <- groups

	for (g in groups) {
		idx = which(groups==g)
		temp <- density(as.numeric(data[data[,DEFAULTS$grouping]==g,'alpha_power']), n = nr_int)
		#temp <- as.numeric(data[data[,DEFAULTS$grouping]==g,'alpha_power'])
		power[idx,] <- as.matrix(temp$y,nrow=1)
	}
	# adjust data format for ggplot
	power <- melt(power)
	# plot density estimation for alpha power distribution
	p.dens <- ggplot(data = power, aes(x=Var2, y= value, linetype=Var1, fill =Var1)) +
		geom_line() + 
		geom_area(alpha=.5, position = 'identity')+
		DEFAULTS$theme +
		scale_fill_manual(values=scoring$colours) +
		theme(legend.position = 'right',
			  axis.title.x=element_text()) +
		labs( title = 'despotism distribution', 
			  y='alpha power',
			  x='estimation intervals')
	# boxplot for alpha power within groups
	p.box <- ggplot(data = power, aes( y= value, fill = Var1)) +
			 geom_boxplot() +
			 DEFAULTS$theme +
			 theme(legend.position = 'right',
			  axis.ticks.x=element_blank(),
			  axis.text.x=element_blank(), 
			  axis.title.x = element_text()) +
			 scale_fill_manual(values=scoring$colours) +
			 labs( title = 'despotism comparison', 
			  y='',
			  x=as.character(DEFAULTS$grouping))
 
 	p <- grid.arrange(p.dens, p.box, ncol=2, widths = 2:1)
 	log_msg(LOG_FILE,paste('PLOTTING:','despotism alpha power', sep='   '))
}

#################
# random forest classification plotting
#################

plot_accuracies <- function(data) {
	acc = as.data.frame(data$accuracies)
	p <- ggplot(data=acc,aes(y=V1)) +
		geom_boxplot(fill=cm.colors(1)) +
		DEFAULTS$theme + 
		theme(axis.text.x=element_blank(), 
			  axis.title.x = element_text()) +
		labs(title='Random Forest accuracies',
			y = 'accuracy')
	return(p)
}

order_data <- function(matrix) {
	mean_order <- colMeans(t(matrix))
	mean_order <- sort(mean_order, decreasing = TRUE)
	data_new <- matrix[match(names(mean_order),rownames(matrix)),]
	return(data_new)
}

get_quant <- function(matrix, perc = c(.9,.75,.5)) {
	quantiles <- matrix(0,nrow=length(perc),ncol=1)
	for (p in perc) {
		quantiles[which(perc==p)] = quantile(matrix,p)
	}
	rownames(quantiles) <- as.character(perc*100)
	return(quantiles)
}


plot_gini <- function(data) {
	data <- order_data(data$gini)
	quantiles <- get_quant(data)
	gini <- melt(data)
	nr_var <- nrow(data)
	p <- ggplot(data=gini) +
		 geom_boxplot(aes(x=Var1, y = value, fill = Var1)) +
		 scale_fill_manual(values=as.matrix(cm.colors(nr_var))) +
		 geom_segment(aes(x=0,xend=nr_var-3,y=quantiles[1],yend=quantiles[1]), linetype = 'dotted') +
		 geom_text(x=nr_var-1, y = quantiles[1], label=paste(rownames(quantiles)[1],'% quantile'), size = 4) +
		 geom_segment(aes(x=0,xend=nr_var-3,y=quantiles[2],yend=quantiles[2]), linetype = 'dotted') +
		 geom_text(x=nr_var-1, y = quantiles[2], label=paste(rownames(quantiles)[2],'% quantile'), size = 4) +
		 geom_segment(aes(x=0,xend=nr_var-3,y=quantiles[3],yend=quantiles[3]), linetype = 'dotted') +
		 geom_text(x=nr_var-1, y = quantiles[3], label=paste(rownames(quantiles)[3],'% quantile'), size = 4) +
		 DEFAULTS$theme +
		 theme(axis.text.x=element_text(angle=90)) +
		 labs(title='Gini Index over classification runs',
		 	  y='gini index',
		 	  x='behaviours')
	return(p)
}

plot_classification <- function(data) {
	p.acc <- plot_accuracies(data)
	p.gini <- plot_gini(data)
	p <- grid.arrange(p.acc,p.gini, ncol=2, widths=1:2)
	log_msg(LOG_FILE,paste('PLOTTING:','random forest results', sep='   '))
	return(p)
}





# batches ordering:
# B5  B2  B6  B10 B7  
# B1  B8  B9  B11 B12



#############################################################
# wrapper for PLOTTING OF ALL ANALYSIS STEPS				#
#############################################################

plotting_all <- function(results, ana_steps) {
	gt_batch <- get_gt_batch_matching(scoring$data)
	if ('sna' %in% ana_steps) {
		# plott sna graphs

		# plot sna parameters
		for (beh in scoring$behaviours) {
			if (nrow(results$sna_parameter[[beh]]) > 0) {
				pdf(paste(PATH$folder[1],'/SNA_parameter_',beh,'.pdf',sep=''), width = 30, height = 15)
					circle_wrapper(results$sna_parameter[[beh]],beh)
				dev.off()
				log_msg(LOG_FILE,paste('PLOTTING:','SNA network parameter', sep='   '))
			}
		}
	}

	if ('glicko' %in% ana_steps) {
		#plot glicko history
		plots <- list()
		for (beh in scoring$behaviours) {
			plots[[beh]] <- glicko_ggplot_wrapper(results$glicko[[beh]],scoring)
			title=paste('A: glicko ranking history for ',beh,sep='')
			xaxistitle='Agonistic interactions'
			yaxistitle=''
			pdf(paste(PATH$folder[1],'/Glicko_history_',beh,'.pdf',sep=''), width = 12, height = 8)
				arange_plots(plots[[beh]],title,yaxistitle,xaxistitle)
			dev.off()
			log_msg(LOG_FILE,paste('PLOTTING:','Glicko ranking history', sep='   '))
		}
		#plot glicko rank correlations
		plots <- list()
		for (beh in scoring$behaviours) {
			temp <- rank_cor_wrapper(results$glicko[[beh]],0.05,scoring)
			title=paste('A: glicko ranking correlation for ',beh,sep='')
			xaxistitle='Intervals of interactions'
			yaxistitle=''
			pdf(paste(PATH$folder[1],'/Glicko_rank_correlations_',beh,'.pdf',sep=''), width = 12, height = 8)
				arange_plots(temp$plots,title,yaxistitle,xaxistitle)
			dev.off()
			log_msg(LOG_FILE,paste('PLOTTING:','Glicko ranking correlation', sep='   '))
		}
	}

	if ('classification' %in% ana_steps) {
		#plot classification accuray over runs and boxplot of gini indices over runs
		pdf(paste(PATH$folder[1], '/Classification_results.pdf', sep = ''), width = 12, height = 6)
			plot_classification(results$classification)
		dev.off()

	}

	if ('despotism' %in% ana_steps) {
		# plot distribution of despotic power via density function over grouping variable
		pdf(paste(PATH$folder[1], '/Despotism_distribution.pdf', sep = ''), width = 12, height = 6)
			plot_despotism(results$despotism)
		dev.off()

	}
}


#############################################################
# wrapper for RANDOM FOREST CLASSIFICATION 					#
#############################################################
run_classification <- function(results, multiple=TRUE) {
	log_msg(LOG_FILE,'RANDOM FOREST: start classification')
	data <- results$frequencies
	parameter <- initialize_param(data)
	log_msg(LOG_FILE,paste('RANDOM FOREST: number of runs: ', as.character(parameter$n.runs)))
	if (multiple) {
		pb <- tkProgressBar(title = "Random Forest computation for multiple runs", min = 0,
                max = parameter$n.runs, width = 300)
		runs <- initialize_multiple_runs(parameter)
		for (r in 1:parameter$n.runs) {
			prediction <- run_cross_validation(data,parameter)
			runs$accuracies[r] <- prediction$accuracy
			runs$gini <- cbind(runs$gini, prediction$gini)
			setTkProgressBar(pb, r, label=paste( round((r/parameter$n.runs)*100), "% completed"))
		}
		close(pb)
		return(runs)
	} else {
		prediction <- run_cross_validation(data,parameter)
		return(prediction)
	}

}
#############################################################
# helper functions for RANDOM FOREST CLASSIFICATION 		#
#############################################################

run_cross_validation <- function(data, parameter) {
	# initialize empty prediction labels
	pt = matrix(,ncol = 1,nrow = parameter$n.animals)
	gini <- matrix(,ncol=parameter$n.animals,nrow = parameter$n.var)
	for ( id in 1:parameter$n.animals) {
		test <- data[id, parameter$variables]
		train <- balance_classes(data[-id,], parameter, id)
		rf <- randomForest(train$data.balanced,as.factor(train$labels.balanced),importance=TRUE)
		pt[id] = rf$classes[predict(rf,test,type = 'response')]
		gini[,id] <- rf$importance[,4]
	}
	rownames(gini) <- parameter$variables
	colnames(gini) <- rownames(data)
	prediction <- list()
	prediction[['labels']] <- pt
	prediction[['gini']] <- gini
	prediction[['accuracy']] <- get_accuracy(prediction,parameter$labels)
	return(prediction)
}

get_accuracy <- function(prediction,groundtruth) {
	classes = unique(groundtruth)
	class.sizes <- table(groundtruth)
	true.labels <- list()
	for (c in classes) {
		true.labels[[c]] <- 0
	}
	for (i in 1:length(prediction$labels)) {
		if (prediction$labels[i] == groundtruth[i]) {
			# convert to within labels instead of count of animnals
			# if labels[i] -< true.labels_i ++1
			if (i <= class.sizes[1]) {
				true.labels[[classes[[1]]]] <- true.labels[[classes[[1]]]] + 1
			} else {
				true.labels[[classes[[2]]]] <- true.labels[[classes[[2]]]] + 1
			}
		}
	}
	acc = ( true.labels[[1]] / class.sizes[[1]] + true.labels[[2]] / class.sizes[[2]] ) / 2
	return(acc)
}


balance_classes <- function(data,parameter,id) {
	train <- list()
	unbalanced_labels = parameter$labels[-id]
	data = data[,parameter$variables]
	class.sizes <- table(unbalanced_labels)
	class.min.id = which(class.sizes == min(class.sizes))
	class.difference = max(class.sizes) - min(class.sizes)
	class.oversample = sample(which(unbalanced_labels == names(class.min.id)))[1:class.difference]
	train[['data.balanced']] <- rbind(data,data[class.oversample,])
	train[['labels.balanced']] <- c(unbalanced_labels,unbalanced_labels[class.oversample])
	return(train)
}

initialize_param <- function(data){
	parameter <- list()
	parameter[['n.animals']] <- nrow(data)
	parameter[['variables']] <- names(data)[-length(names(data))]
	parameter[['n.var']] <- length(parameter$variables)
	parameter[['n.runs']] <- 10
	parameter[['labels']] <- data$genotype
	return(parameter)
}

initialize_multiple_runs <-function(parameter) {
	runs <- list()
	runs[['accuracies']] <- matrix(,nrow = parameter$n.runs, ncol=1)
	runs[['gini']] <- matrix(,ncol=0,nrow = parameter$n.var)
	rownames(runs$gini) <- parameter$variables
	return(runs)
}

########
#
# ADD:    possibility to add variables or use different data
#############################################################
# wrapper for additional metrics and classification input	#
#############################################################

# count for local place preference of animals
get_place_preference <- function(scoring) {
	log_msg(LOG_FILE, paste('PLACE PREFERENCE: computing place preference'))
	places <- unique(scoring$mapped.data$place)
	batches <- scoring$batches
	temp <- as.data.frame(matrix(,ncol = length(places), nrow = 0))
	colnames(temp) <- places
	for (b in batches) {
		preference <- table(scoring$mapped.data[scoring$mapped.data$batch == b,]$place)
		temp <- bind_rows(temp,preference)
	}
	rownames(temp) <- batches
	temp[is.na(temp)] <- 0
	return(temp)
}

## get behaviour frequencies per animal as classification input
get_frequencies <- function(scoring, phase_split=TRUE) {
	log_msg(LOG_FILE, paste('BEHAVIOUR FREQUENCIES: computing behaviour occurences'))
	scoring$data <- remove_empty_behaviour(scoring$data)
	animals <- unique(scoring$data[[DEFAULTS$identifier]])
	behaviour <- unique(scoring$data$behaviour)
	freq <- as.data.frame(matrix(,ncol = length(behaviour)+1, nrow = 0))
	colnames(freq) <- c(behaviour,'genotype')
	#####
	# add split by phase
	for (i in animals) {
		temp = table(scoring$data[scoring$data[[DEFAULTS$identifier]] == i,]$behaviour)
		temp[['genotype']] = scoring$data[scoring$data[[DEFAULTS$identifier]] == i,]$genotype[1]
		freq = bind_rows(freq,temp)
	}
	rownames(freq) <- animals
	freq[is.na(freq)] <- 0
	colnames(freq) <- c(behaviour,'genotype')
	return(freq)
}

#############################################################
# DESPOTISM ANALYSIS TO INCLUDE!!!!!!!!!!					#
#############################################################


get_despotism <- function(glicko) {
	nr_updates <- dim(glicko$history)[2]
	nr_animals <- dim(glicko$history)[1]
	final_rating = sort(glicko$history[,nr_updates,'Rating'])
	alpha_beta_power = final_rating[nr_animals]-final_rating[nr_animals-1]
	alpha_max_power = final_rating[nr_animals]-final_rating[1]
	imposed_power = alpha_beta_power / alpha_max_power
	return(imposed_power)
}

despotism_all <- function(glicko,scoring) {
	gt_batch <- get_gt_batch_matching(scoring$data)
	despotism <- matrix(,nrow = length(gt_batch$batch), ncol=2)
	colnames(despotism) <- c('alpha_power',as.character(DEFAULTS$grouping))
	rownames(despotism) <- gt_batch$batch
	for (b in gt_batch$batch) {
		log_msg(LOG_FILE,paste('DESPOTISM: computing despotic power for batch',as.character(b)))
		temp = glicko[[1]][[b]]
		despotism[which(gt_batch$batch == b),1] = get_despotism(temp)
		despotism[which(gt_batch$batch == b),2] = as.character(gt_batch$genotype[gt_batch$batch==b])
		}
	return(despotism)
	# write.csv(as.data.frame(despotism), paste(Sys.Date(),'_despotism_ranking.csv',sep=''), row.names = TRUE)
}

#############################################################
# wrapper for input and analysis steps						#
#############################################################


get_input <- function() {
	win.select <- tktoplevel()
	tktitle(win.select) = 'Select working directory and analysis steps'

	cb1 <- tkcheckbutton(win.select) # SNA check box
	cb2 <- tkcheckbutton(win.select) # glicko check box
	cb3 <- tkcheckbutton(win.select) # random forest check box
	cb4 <- tkcheckbutton(win.select) # despotism check box
	cb5 <- tkcheckbutton(win.select) # plotting of results
	# initial values unselect all steps
	cb1Value <- tclVar("0")
	cb2Value <- tclVar("0")
	cb3Value <- tclVar("0")
	cb4Value <- tclVar("0")
	cb5Value <- tclVar("0")
	tkconfigure(cb1,variable=cb1Value)
	tkconfigure(cb2,variable=cb2Value)
	tkconfigure(cb3,variable=cb3Value)
	tkconfigure(cb4,variable=cb4Value)
	tkconfigure(cb5,variable=cb5Value)
	tkgrid(tklabel(win.select,text="social network analysis"),cb1)
	tkgrid(tklabel(win.select,text="glicko rating analysis"),cb2)
	tkgrid(tklabel(win.select,text="random forest classification"),cb3)
	tkgrid(tklabel(win.select,text="despotism analysis"),cb4)
	tkgrid(tklabel(win.select,text="plotting of results"),cb5)
	Okquit <- function() {
		if (sum(as.numeric(tclvalue(cb1Value)),as.numeric(tclvalue(cb2Value)), as.numeric(tclvalue(cb3Value)), as.numeric(tclvalue(cb4Value))) == 0) {
		tkmessageBox(message= 'Please select steps to perform')
		} else { tkdestroy(win.select)}

	}
	OK.but <- tkbutton(win.select,text="OK",command=Okquit)
	tkgrid(OK.but)
	tkfocus(win.select)
	# stop program until seperator is selected
	tkwait.window(win.select)
	# assign selected values to return variable
	run = list()
	run$sna = as.numeric(tclvalue(cb1Value))
	run$glicko = as.numeric(tclvalue(cb2Value))
	run$classification = as.numeric(tclvalue(cb3Value))
	run$despotism = as.numeric(tclvalue(cb4Value))
	run$plotting = as.numeric(tclvalue(cb5Value))
	return(run)
}

# initialize result files for selected analysis steps
init_results <- function(run) {
	results <- list()
	for (ana in names(run)) {
		if (run[[ana]] == 1) {
			results[[ana]] <- list()
		}
	}
	return(results)
}

run_analysis <- function(scoring, run, export = FALSE) {
	#initialize result list
	results <- init_results(run)
	# map beehaviours to apriori deefined groups
	scoring <- behave_mapping(scoring)
	scoring$batches <- get_batches(scoring$data)
	# data$data <- color_mapping(data$data,batches)
	scoring <- check_encoding(scoring)
	animals <- get_animals(scoring)
	gt_batch <- get_gt_batch_matching(scoring$data)
	# add additional metrics
	results[['place_preference']] <- get_place_preference(scoring)
	results[['frequencies']] <- get_frequencies(scoring)
	# get steps to perform
	ana_steps <- names(run)[run==1]
	log_msg(LOG_FILE, paste('SELECTED ANALYSIS STEPS: ',paste(ana_steps, collapse = ' - '),sep=''))
	# run selected analysis steps
	for (beh in scoring$behaviours) {
		temp.behav <- behav_split(scoring$mapped.data, beh)

		for (b in scoring$batches) {
			temp <- batch_split(temp.behav, b)
			
			if (run$sna) {
				if (nrow(temp) == 0) {
					results$sna[[beh]][[b]] <- 'no given interaction'
				} else {
					results$sna[[beh]][[b]] <- create_networks(data_split(temp, DEFAULTS$split), animals[b])
					if (export) {
						SNA <- results$sna
						save(SNA, file=paste(PATH$folder[2],'/Behavioural_networks_', Sys.Date(),'.RData', sep = ''))
					}
				}
			log_msg(LOG_FILE, paste('SNA: computing behaviour',as.character(beh),'for batch',as.character(b)))
			}
			if (run$glicko) {
				if (nrow(temp) == 0) {
					results$glicko[[beh]][[b]] <- 'no given interaction'
				} else {
					results$glicko[[beh]][[b]] <- check_interactions(temp)
					if (export) {
						GLICKO <- results$glicko
						save(GLICKO, file=paste(PATH$folder[2],'/Behavioural_glicko_', Sys.Date(),'.RData', sep = ''))
					}
				}
			# include logging
			log_msg(LOG_FILE, paste('GLICKO: computing behaviour',as.character(beh),'for batch',as.character(b)))
			}
		}
		if (run$sna) {
			results$sna_parameter[[beh]] <- get_parameter(results$sna[[beh]], gt_batch)
		log_msg(LOG_FILE, paste('SNA-TOPOLOGY: computing network parameter for ',as.character(beh)))
			if (export) {
				SNA.PARAM <- results$sna_parameter
				save(SNA.PARAM, file=paste(PATH$folder[2],'/Behavioural_network_parameter', Sys.Date(),'.RData', sep = ''))
			}
		}
	}
		# run classification based on behaviour frequencies
	if (run$classification) {
		results[['classification']] <- run_classification(results)
		if (export) {
			CLASSIFICATION <- results$classification
			save(CLASSIFICATION, file=paste(PATH$folder[2],'/Behavioural_classification_', Sys.Date(),'.RData', sep = ''))
			}	
	}
	if (run$despotism) {
		results[['despotism']] <- despotism_all(results$glicko, scoring)
		if (export) {
			DESPOTISM <- results$despotism
			save(DESPOTISM, file=paste(PATH$folder[2],'/Despotism_analysis_', Sys.Date(),'.RData', sep = ''))
			}	
	}
	if (run$plotting) {
		plotting_all(results, ana_steps)
	}
	if (export) {
		PLACEPREFERENCE <- results$place_preference
		save(PLACEPREFERENCE, file=paste(PATH$folder[2],'/Behavioural_place_preference_', Sys.Date(),'.RData', sep = ''))
		FREQUENCIES <- results$frequencies
		save(FREQUENCIES, file=paste(PATH$folder[2],'/Behaviour_frequencies_', Sys.Date(),'.RData', sep = ''))
	}

	#####

	return(results)
}

