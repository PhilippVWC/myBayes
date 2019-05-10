### global variables and constants ###
### P.v.W.Crommelin ###
N_ALPHA = 100 # number of grid ticks for a
N_R = 200 # number of grid ticks for r
N_PRED = 200 # number of iterations of time series
N_DATA = 20 # number of iterations of time series
SIGMA = 0.01
X0 = 0.3 # starting value
R = 0.1
ALPHA = 2.0

# Graphical parameters
DASHDOTTED_LONG = 6
DASHED_LONG = 5
DASHDOTTED = 4
DOTTED = POINTS = 3 # Linetype for plots
DASHED = 2
SOLID = 1
LWD_GRID = 0.75 #line width for grids
COLOR_BG_LEG = "#FFFF0099" # transparent yellow
COLOR_BG = "#00005599"
COLOR_AXIS = colors()[41] # light brown
WIDTH = 1000
HEIGHT = 400
NX_GRID = 10
NY_GRID = 10
COLOR_GRID = "grey"
PARALLEL = 0
HORIZONTAL = 1
PERPENDICULAR = 2


###### enable certain plots ######
PLOT_MAP = TRUE
SAVE_MAP = TRUE
PLOT_LIKELIHOOD = FALSE
SAVE_LIKELIHOOD = TRUE
PLOT_TIME_SERIES = TRUE
SAVE_TIME_SERIES = TRUE
PLOT_STD = TRUE
SAVE_STD = TRUE
PLOT_FORECAST = TRUE
SAVE_FORECAST = TRUE