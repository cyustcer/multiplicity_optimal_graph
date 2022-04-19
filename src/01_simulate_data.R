.libPaths("/CHBS/apps/EB/software/R-bundle-Novartis/0.1-gomkl-2019a-R-4.1.0")
library(optparse, lib.loc='/home/chenya3i/Packages/R')
library(here)
source(here::here("src/utils/grids.R"))
source(here::here("src/utils/simulation.R"))

# Parse simulation parameters
opt_list = list(make_option("--number", type="integer", default=1),
                make_option("--correlation", type="double", default=0.0),
                make_option("--n_simulation", type="integer", default=10^6),
                make_option("--weight", type="character", default="1/3, 1/3, 1/3"),
                make_option("--type_1_error", type="double", default=0.025),
                make_option("--margin_power", type="character", default="0.8, 0.8, 0.8"),
                make_option("--w_step", type="double", default=0.1),
                make_option("--G_step", type="double", default=0.1),
                make_option("--w_range", type="character", default="0, 1"),
                make_option("--batch_size", type="integer", default=1000))

opt_parser = OptionParser(option_list=opt_list)
opt = parse_args(opt_parser)


# Convert parameters
n = opt$number
corr = c(opt$correlation)
n.sim = opt$n_simulation
eval(parse(text=paste0("weight = c(", opt$weight, ")")))
tie = opt$type_1_error
eval(parse(text=paste0("mp = c(", opt$margin_power, ")")))
n.hypo = length(mp)
sigma = matrix(corr, nrow=n.hypo, ncol=n.hypo)
diag(sigma) = 1
w.step = opt$w_step
G.step = opt$G_step
eval(parse(text=paste0("w.range = c(", opt$w_range, ")")))
bs = opt$batch_size


# Create grid
data.grid <- gen.data.grid(n.hypo, w.step=w.step, w.range=w.range, G.step=G.step)
grids <- load.grid(data.grid)
ws = grids$ws
Gs = grids$Gs


# Extract grid batch
ws = ws[((n-1)*bs+1):min(n*bs, dim(data.grid)[1]), ]
Gs = Gs[, , ((n-1)*bs+1):min(n*bs, dim(data.grid)[1])]


# Data simulation and save
data.sim <- simulate.grid.power(n.sim, mp, ws, Gs, sigma, tie, weight)
file.dir = here::here("outputs/grids", paste0("mp-", paste(mp, collapse="-"),
                                              "_weight-", paste(weight, collapse="-"),
                                              "_corr-", corr, "_n.sim-", log(n.sim)/log(10),
                                              "_w.step-", w.step, "_w.range-", paste(w.range, collapse="-"),
                                              "_G.step-", G.step))
file.name = paste0("grid.index-", n, ".RData")
if (!dir.exists(file.dir)) {
  dir.create(file.dir)
}
save(data.sim, file=here::here(file.dir, file.name))