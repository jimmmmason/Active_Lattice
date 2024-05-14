#Set up
cd("/home/jm2386/Active_Lattice/");
# import Pkg; Pkg.add("DrWatson")
using DrWatson;
@quickactivate "Active_Lattice";
include("/home/jm2386/Active_Lattice/src/pm_pdes.jl");
include("/home/jm2386/Active_Lattice/src/pm_plot.jl");

# Access the command-line argument 'i'
try
    input = parse(Int64, ARGS[1]);
    i_max = 3
    i = (input %i_max) +1
    j = (input ÷ i_max)+1
    println("$i, $j")
catch
    println("no input")
end
# wait 1min
sleep(60)

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=some.email@your_institutions.com