#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))

usage = '\nPermutation test to determine whether two groups are significantly
different across longitudinal data that may have differing patterns
of data for each subject (time, n, etc). Applies a permutation test using
splines fit to the two groups. Categorical labels are shuffled among the
subjects to determine a random distribution for comparison.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table, prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest; must be a column header in table',
              default=NA, type = 'character'),
  make_option(c('--groups'),
              help='If >2 groups in category, the two of interest: comma-delimited (e.g. groups=Healthy,Sick)',
              default=NA, type = 'character'),
  make_option(c('-x', '--x_variable'),
              help='REQUIRED: The independent variable (e.g. time); must be column header',
              default=NA, type = 'character'),
  make_option(c('-y', '--y_variable'),
              help='REQUIRED: The response variable; must be a column header',
              default=NA, type = 'character'),
  make_option(c('-p', '--unit_id'),
              help='REQUIRED: The column header for your grouping (e.g. Patient_ID, User_Name, etc)',
              default=NA, type = 'character'),
  make_option(c('--perms'),
              help='Number of permutation shuffles [default %default]',
              default=999, type = 'integer'),
  make_option(c('--cut'),
              help='Cut Unit IDs that occur less than this many times',
              default=NA, type = 'character'),
  make_option(c('--intervals'),
              help='Number of sampling intervals along spline [default %default]',
              default=10000, type = 'integer'),
  make_option(c('--spar'),
              help='The spar parameter when fitting splines (0 - 1); default is calculated from data',
              default=NULL),
  make_option(c('--plot'),
              help='Plot the data too! Provide a filename/path, with .png extension',
              default=NA, type = 'character')
)
opt = parse_args(OptionParser(usage=usage, option_list=option_list))

if (is.na(opt$input) | is.na(opt$category) | is.na(opt$x_variable) |
  is.na(opt$y_variable) | is.na(opt$unit_id)) {
    stop('Missing required parameters. See usage and options (--help)')
  }

# args = commandArgs(trailingOnly=TRUE)

# Parse command line
infile = opt$input  # tsv file with data in long form
category = opt$category  # the column header label for the two groups
groups = opt$groups  # the two groups of interest, if column has >2
x.cat = opt$x_variable  # the time series label
y.cat = opt$y_variable  # the response variable label
unit.id = opt$unit_id  # the header label defining the individuals (patients, etc)
num.perm = as.numeric(opt$perms)  # default 999
cut.low = opt$cut
spar.param = opt$spar # default NULL
samp.intervals = opt$intervals 
plot.results = opt$plot  # name of the plot file.png
shuff.id = 'cat.shuff'

# Read infile
df = read.delim(file = infile, header = 1, check.names = F, sep = '\t')
if (is.na(groups)) {
  if (length(unique(df[, category])) > 2) {
    stop('More than two groups in category column. Define groups with "--groups=Name1,Name2"')
  }
  v1 = unique(df[, category])[1]
  v2 = unique(df[, category])[2]
} else {
  v1 = strsplit(groups, ',')[[1]][1]
  v2 = strsplit(groups, ',')[[1]][2]
  }

if (!is.na(cut.low)) {
  cut_low = as.numeric(cut.low)
  keep.ids = data.frame(table(df[, unit.id]))
  keep.ids = as.character(keep.ids[keep.ids$Freq > cut.low, ]$Var1)
  df = df[df[,unit.id] %in% keep.ids, ]
}


cat(paste('\nTesting between', v1, 'and', v2, 'for a difference in', y.cat, '\n'))
cat(paste('\nScalpel please: performing permusplinectomy with', num.perm, 'permutations...\n'))

# The experimentally reported response
df.v1 = df %>% filter(df[, category] == v1 & !is.na(df[, x.cat]))
df.v2 = df %>% filter(df[, category] == v2 & !is.na(df[, x.cat]))
df.v1.spl = with(df.v1,
                smooth.spline(x=df.v1[, x.cat], y=df.v1[, y.cat],
                spar = spar.param))
df.v2.spl = with(df.v2,
                smooth.spline(x=df.v2[, x.cat], y=df.v2[, y.cat],
                spar = spar.param))
x0 = max(c(min(df.v1.spl$x)), min(df.v2.spl$x))
x1 = min(c(max(df.v1.spl$x)), max(df.v2.spl$x))
xby = (x1 - x0) / samp.intervals
xx = seq(x0, x1, by = xby)
v1.spl.f = data.frame(predict(df.v1.spl, xx))
colnames(v1.spl.f) = c('x', 'var1')
v2.spl.f = data.frame(predict(df.v2.spl, xx))
colnames(v2.spl.f) = c('x', 'var2')
real.spl.dist = merge(v1.spl.f, v2.spl.f, by = 'x')
real.spl.dist$abs.distance = abs(real.spl.dist$var1 - real.spl.dist$var2)
real.area = sum(real.spl.dist$abs.distance) / samp.intervals

# Define the permutation function
spline_permute = function(randy, unit.id, category, x.cat, y.cat) {
  randy.meta = randy %>% distinct_(unit.id, .keep_all = T)
  randy.meta$cat.shuff = sample(randy.meta[,category])
  randy.meta = randy.meta %>% select_(unit.id, shuff.id)
  randy = merge(randy, randy.meta, by = unit.id, all = T)
  randy.v1 = filter(randy, cat.shuff == v1 & !is.na(randy[,x.cat]))
  randy.v2 = filter(randy, cat.shuff == v2 & !is.na(randy[,x.cat]))
  randy.v1.spl = with(randy.v1,
                      smooth.spline(x=randy.v1[, x.cat], y=randy.v1[, y.cat]))
  randy.v2.spl = with(randy.v2,
                      smooth.spline(x=randy.v2[, x.cat], y=randy.v2[, y.cat]))
  x0 = max(c(min(randy.v1.spl$x)), min(randy.v2.spl$x))
  x1 = min(c(max(randy.v1.spl$x)), max(randy.v2.spl$x))
  xby = (x1 - x0) / samp.intervals
  xx = seq(x0, x1, by = xby)
  randy.v1.fit = data.frame(predict(randy.v1.spl, xx))
  colnames(randy.v1.fit) = c('x', 'var1.y')
  randy.v2.fit = data.frame(predict(randy.v2.spl, xx))
  colnames(randy.v2.fit) = c('x', 'var2.y')
  spl.dist = merge(randy.v1.fit, randy.v2.fit, by = 'x')
  spl.dist$abs.distance = abs(spl.dist$var1.y - spl.dist$var2.y)
  perm.area = sum(spl.dist$abs.distance) / samp.intervals
  permuted = append(permuted, perm.area)
  return(permuted)
}

# Run the permutation over desired number of iterations
permuted = list()
permuted = replicate(num.perm, 
                     spline_permute(df, unit.id, category, x.cat, y.cat))
pval = (sum(permuted >= as.numeric(real.area)) + 1) / (num.perm + 1)

# Return the p-value
cat(paste('\np-value =', round(pval, digits = 5), '\n\n'))
spar.param = NULL
if (!is.na(plot.results)) {
  df.p = rbind(df.v1, df.v2)
  df.pick = c(x.cat, category, y.cat)
  plot.df = df.p[, df.pick]
  plot.df = plot.df[!is.na(plot.df[, x.cat]), ]
  plot.df = droplevels(plot.df)
  p = ggplot(plot.df, aes(x=plot.df[,x.cat], y=plot.df[,y.cat], color=as.character(plot.df[,category]))) +
    geom_point() + geom_smooth(span = spar.param) + xlab(x.cat) + ylab(y.cat) +
    scale_color_manual(name=category, values = c("#0072B2","#D55E00"))
  # p
  ggsave(plot.results,
         height = 3.5, width = 4, units = 'in', dpi = 600)
  cat(paste('Plot saved\n'))
}


