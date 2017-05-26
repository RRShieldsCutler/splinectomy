#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(reshape2))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(cowplot))

# TODO:
  # Optparsing
  # Remove cowplot dependency
  # Remove dplyr dependency

usage = '\nSliding spline test to compare a binary categorical variable across
  longitudinal data that may be sparse or messy. Generates splines for each
  individual and uses those points to compare groups across time. Result is
  a set of non-parametric p-values at user-defined density across the time
  series. Returns spline plot, p-value plot (png), and p-value table (txt) to
  the current working directory.'

option_list = list(
  make_option(c('-i', '--input'),
              help='REQUIRED: Input is long-form tab-delimited table, prefiltered for missing values',
              default=NA, type = 'character'),
  make_option(c('-c', '--category'),
              help='REQUIRED: Categorical variable of interest; must be a column header in table',
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
  make_option(c('--spar'),
              help='The spar parameter when fitting splines (0 - 1); default is calculated from data',
              default=NULL),
  make_option(c('--spline_intervals'),
              help='Number of intervals extrapolated across the splines [default %default]',
              default=100, type = 'integer'),
  make_option(c('--density'),
              help='Data points required for each sliding spline position to test [default %default]',
              default=3, type = 'integer'),
  make_option(c('--n_per_unit'),
              help='Data points required per unit (e.g. patient) to keep in dataset [default %default]',
              default=3, type = 'integer'),
  make_option(c('--prefix'),
              help='Prefix for results files',
              default='', type = 'character')
)
opt = parse_args(OptionParser(usage=usage, option_list=option_list))

# Parse command line
infile = opt$input  # tsv file with data in long form
category = as.character(opt$category)  # the column header label for the two groups
x.cat = opt$x_variable  # the time series label
y.cat = opt$y_variable  # the response variable label
unit.id = opt$unit_id  # the header label defining the individuals (patients, etc)
num.bits = as.numeric(opt$spline_points)  # default 100
sparsity = as.numeric(opt$density)  # cutoff for number of points at given time; default 3 (>= 3)
unit_number = as.numeric(opt$n_per_unit)  # min number of points needed per individual
spar.param = opt$spar  # default NULL
prefix = opt$prefix

## DEBUGGING TEST PARAMETERS
# setwd('~/Box Sync/knights_box/bile_vsg/data/')
# infile = 'pre_surg_vitals_relweight_response.txt'  # tsv file with data in long form
# category = as.character('response')  # the column header label for the two groups
# x.cat = 'YEARS_REL_TO_SURG'  # the time series label
# y.cat = 'WEIGHT_REL_TO_PRE_SURG'  # the response variable label
# unit.id = 'PATIENT_ID'  # the header label defining the individuals (patients, etc)
# num.bits = as.numeric(100)  # default 100
# sparsity = as.numeric(4)

cat(paste('Running sliding spline test with', num.bits,
          'time points extrapolated from splines...\n'))

# Read infile and limit to ids matching sparsity parameter
df = read.delim(file = infile, header = 1, check.names = F, sep = '\t')
unit.ids.tab = data.frame(table(df[, unit.id]))
unit.ids.notsparse = unit.ids.tab[unit.ids.tab$Freq >= unit_number, ]
unit.ids.keep = as.character(unit.ids.notsparse$Var1)
df = df[df[, unit.id] %in% unit.ids.keep, ]

v1 = unique(df[, category])[1]
v2 = unique(df[, category])[2]

# Get the range of the independent variable (e.g. time) to set spline limits
x.min = min(df[, x.cat])
x.max = max(df[, x.cat])
xx = seq(x.min, x.max, by = ((x.max - x.min) / (num.bits - 1)))
spl.table = setNames(data.frame(xx), c('x'))

# Save the group labels for each individual/unit
df.groups = df %>% distinct_(unit.id, .keep_all = T)
df.groups = df.groups %>% select_(unit.id, category)

# Generate splines for each individual
for (i in unit.ids.keep) {
  unit.df = subset(df, df[, unit.id]==i)
  unit.spl = with(unit.df,
                  smooth.spline(x=unit.df[, x.cat], y=unit.df[, y.cat],
                  spar = spar.param))
  xx.i = subset(xx, xx >= min(unit.spl$x) & xx <= max(unit.spl$x))
  unit.spl.f = data.frame(predict(unit.spl, xx.i))
  colnames(unit.spl.f) = c('x', i)
  spl.table = merge(spl.table, unit.spl.f, by = 'x', all = T)
}

# Prepare the spline table for statistical testing
spl.table.p = tibble::column_to_rownames(df = spl.table, var = 'x')
spl.table.p = as.data.frame(t(spl.table.p))
spl.table.p = tibble::rownames_to_column(spl.table.p, var = unit.id)
spl.table.p = merge(spl.table.p, df.groups, by = unit.id, all = T)

# Define the non-parametric test function
mann_whitney_per_bit = function(spline.table, n) {
  n = as.character(n)
  x.pick = c(unit.id, n, category)
  x.table = spline.table[, x.pick]
  x.table = x.table[!is.na(x.table[, n]), ]
  cat.freq = data.frame(table(x.table[, category]))
  if (cat.freq$Freq[1] >= sparsity & cat.freq$Freq[2] >= sparsity) {
    mw = wilcox.test(x.table[, n] ~ x.table[, category])
    pval = mw$p.value
    num.pts = cat.freq$Freq[1] + cat.freq$Freq[2]
    return(c(pval, n, num.pts))
  } else {}
}

# Run the stats test on each step of the spline on each data unit
pvals.list = list()
pval.x = list()
pval.num.pts = list()  # Save the number of total data points in each comparison
for (n in xx) {
  pval = mann_whitney_per_bit(spl.table.p, n)
  pvals.list = append(pvals.list, pval[1])
  pval.x = append(pval.x, pval[2])
  pval.num.pts = append(pval.num.pts, pval[3])
}
pval.df = do.call(rbind, Map(data.frame,
                x.series=as.numeric(pval.x),
                p.value=as.numeric(pvals.list),
                N=as.numeric(pval.num.pts)))

# Plot the splines
plot.spline.data = melt(data = spl.table, id.vars = 'x')
colnames(plot.spline.data) = c('x', unit.id, 'value')
plot.spline.data = merge(plot.spline.data, df.groups, by = unit.id, all = T)
colnames(plot.spline.data) = c('UNIT', 'x', 'value', 'category')
p = ggplot(plot.spline.data, aes(x=x, y=value, group=UNIT, color=category)) +
  geom_line(na.rm = T) + scale_color_manual(name=category,
  values = c("#0072B2","#D55E00")) +
  xlab(x.cat) + ylab(y.cat) +
  theme(legend.position='right')
# p
ggsave(paste(prefix, 'spline_plot_result.png', sep = ''), height = 3.5, width = 5, units = 'in', dpi = 600)

# Plot the p-values as function of the independent variable (e.g. time)
# Scale the size of the line and points according to the number of observations
norm.range = function(x){(x-min(x)) / (max(x)-min(x))}
pval.df$N.norm = norm.range(pval.df$N)
p = ggplot(pval.df, aes(x=x.series, y=p.value, size=N.norm)) + geom_line() +
  geom_point(shape = 20) +
  geom_hline(aes(yintercept = 0.05), linetype='dashed') +
  xlab(x.cat) + ylab('Mann-Whitney p-value') +
  theme(legend.position='none') + ylim(0,1)
# p
ggsave(paste(prefix, 'spline_evaluated_pvalues.png', sep = ''),
       height = 3.5, width = 4, units = 'in', dpi = 600)

# Save the p-values table to file
pval.df = pval.df[,1:3]
colnames(pval.df) = c(x.cat, 'p-value', 'number of observations')
write.table(pval.df, file = paste(prefix, 'spline_evaluated_pvalues.txt', sep = ''), row.names = F,
            sep = '\t', quote = F)
