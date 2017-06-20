#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(cowplot))

usage = '\nPermutation to test whether there is a non-zero trend among a set
of individuals/samples over a continuous variable (such as time). So, there
does not need to be two groups in this test. The x variable datapoints are
permuated within each case/individual, thus maintaining the distribution in
the y component but shuffling the hypothesized trend.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table, prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest; must be a column header in table',
              default=NA, type = 'character'),
  make_option(c('--group'),
              help='If >1 group in category, the group of interest (e.g. group=disease)',
              default=NA, type = 'character'),
  make_option(c('-x', '--x_variable'),
              help='REQUIRED: The independent variable (e.g. time); must be column header',
              default=NA, type = 'character'),
  make_option(c('-y', '--y_variable'),
              help='REQUIRED: The response variable; must be a column header',
              default=NA, type = 'character'),
  make_option(c('-p', '--unit_id'),
              help='REQUIRED: The column header for your case grouping (e.g. Patient_ID, User_Name, etc)',
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

# Parse command line
infile = opt$input  # tsv file with data in long form
category = opt$category  # the column header label for the group
test_grp = opt$group  # the group of interest, if column has >1
x.cat = opt$x_variable  # the time series label
y.cat = opt$y_variable  # the response variable label
unit.id = opt$unit_id  # the header label defining the individuals (patients, etc)
num.perm = as.numeric(opt$perms)  # default 999
cut.low = opt$cut
spar.param = opt$spar # default NULL
samp.intervals = opt$intervals 
plot.results = opt$plot  # name of the plot file.png
shuff.id = 'y_shuff'


## DEBUGGING DATA
# setwd('~/Box Sync/knights_box/splinectomy/test/')
# infile = '../test/ChickWeight.txt'  # tsv file with data in long form
# category = 'Diet'  # the column header label for the group
# test_grp = '1'  # the group of interest, if column has >1
# x.cat = 'Time'  # the time series label
# y.cat = 'weight'  # the response variable label
# unit.id = 'Chick'  # the header label defining the individuals (patients, etc)
# num.perm = as.numeric(99)  # default 999
# cut.low = NA
# spar.param = NULL # default NULL
# samp.intervals = 100 
# plot.results = 'trendyplot_tests.png'  # name of the plot file.png
# shuff.id = 'y_shuff'


# Read infile
df = read.delim(file = infile, header = 1, check.names = F, sep = '\t')
if (is.na(test_grp)) {
  if (length(unique(df[, category])) > 1) {
    stop('More than one group in category column. Define group with "--group=Name1"')
  }
  v1 = unique(df[, category])[1]
} else {
  v1 = as.character(test_grp)
}
if (!is.na(cut.low)) {
  cut_low = as.numeric(cut.low)
  keep.ids = data.frame(table(df[, unit.id]))
  keep.ids = as.character(keep.ids[keep.ids$Freq > cut.low, ]$Var1)
  df = df[df[,unit.id] %in% keep.ids, ]
}


cat(paste('\nTesting group', v1, 'for a non zero trend in', y.cat, '\n'))
cat(paste('\nPerforming the trendyspline test with', num.perm, 'permutations...\n'))

# The experimentally reported response
df.v1 = df %>% filter(df[, category] == v1 & !is.na(df[, x.cat]))

## Doing the group mean because most cases will likely be asking about the
### behavior of a _group_. i.e. does this group change consistently. If we
### measured individual splines, noisey data could produce false positives

## Then spline for the group
## Measure distance to the mean line
## Within the units, shuffle the y values keeping the set of x
## Sig equals the fraction that are equal or greater than truth

## First determine the group mean (null hypothesis for changing over time)
y.mean = mean(df.v1[, y.cat])
df.v1.spl = with(df.v1,
                 smooth.spline(x=df.v1[, x.cat], y=df.v1[, y.cat],
                               spar = spar.param))
x0 = min(df.v1.spl$x)
x1 = max(df.v1.spl$x)
xby = (x1 - x0) / (samp.intervals - 1)
xx = seq(x0, x1, by = xby)
v1.spl.f = data.frame(predict(df.v1.spl, xx))
colnames(v1.spl.f) = c('x', 'var1')
real.spl.dist = v1.spl.f
real.spl.dist$y_mean = y.mean
real.spl.dist$abs.distance = abs(real.spl.dist$var1 - real.spl.dist$y_mean)
real.area = sum(real.spl.dist$abs.distance) / samp.intervals

# Define the permutation function
spline_permute = function(randy, unit.id, category, x.cat, y.cat) {
  randy.meta = randy %>% select_(unit.id, x.cat, category)
  randy.meta$y_shuff = sample(randy[, y.cat])
  # randy.meta = randy.meta %>% select_(unit.id, shuff.id)
  # randy = merge(randy, randy.meta, by = unit.id, all = T)
  randy.meta = randy.meta %>% filter(!is.na(randy.meta[, x.cat]))
  randy.v1.spl = with(randy.meta,
                      smooth.spline(x=randy.meta[, x.cat], y=randy.meta[, shuff.id]))
  x0 = min(randy.v1.spl$x)
  x1 = max(randy.v1.spl$x)
  xby = (x1 - x0) / (samp.intervals - 1)
  xx = seq(x0, x1, by = xby)
  randy.v1.fit = data.frame(predict(randy.v1.spl, xx))
  colnames(randy.v1.fit) = c('x', 'var1.y')
  spl.dist = randy.v1.fit
  spl.dist$y_mean = y.mean
  spl.dist$abs.distance = abs(spl.dist$var1.y - spl.dist$y_mean)
  perm.area = sum(spl.dist$abs.distance) / samp.intervals
  permuted = append(permuted, perm.area)
  return(permuted)
}

# Run the permutation over desired number of iterations
permuted = list()
permuted = replicate(num.perm, 
                     spline_permute(randy = df.v1, unit.id, category, x.cat, y.cat))
pval = (sum(permuted >= as.numeric(real.area)) + 1) / (num.perm + 1)

# Return the p-value
cat(paste('\np-value =', round(pval, digits = 5), '\n\n'))

if (!is.na(plot.results)) {
  df.p = df.v1
  df.pick = c(x.cat, category, y.cat)
  plot.df = df.p[, df.pick]
  plot.df = plot.df[!is.na(plot.df[, x.cat]), ]
  plot.df = droplevels(plot.df)
  p = ggplot(plot.df, aes(x=plot.df[,x.cat], y=plot.df[,y.cat], color=as.character(plot.df[,category]))) +
    geom_point() + geom_smooth(span = spar.param) + xlab(x.cat) + ylab(y.cat) +
    scale_color_manual(name=category, values = c("#0072B2")) +
    geom_hline(yintercept = y.mean, size = 1.1, color = 'dark red', linetype = 'dashed')
  p
  ggsave(plot.results,
         height = 3.5, width = 4, units = 'in', dpi = 600)
  cat(paste('Plot saved\n'))
}


