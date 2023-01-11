#################################################
# MPS transfer function			#
#################################################
# Example for the paper 
# "On the receive path calibration of Magnetic
#  Particle Imaging Systems" 
# DOI: 10.1109/TIM.2022.3219461
#################################################
# Sample: 10µL perimag 8.5 mg/ml 
# MPS: MPS1 University Medical Center Hamburg-Eppendorf, Hamburg, Germany
# Julia: 1.8

##############################
# Activate local environment #
##############################
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# load packages
using MPIFiles
using Plots, Plots.PlotMeasures
using LazyArtifacts
using Scratch

pyplot() # use PyPlot backend for plotting

squeeze(A) = dropdims(A, dims=tuple(findall(([size(A)...].==1))...))


# Download Data

folderName = artifact"MPSTransferFunction"
const tmpdir  = @get_scratch!("MPSTF")
@info "If you want to check the output of the tests, please head to $tmpdir."
cp(joinpath(folderName,"example.mdf"),joinpath(tmpdir,"example.mdf"),force=true)
chmod(tmpdir,0o777,recursive=true)  # write rights for setting the TF


#####################
# Parameter setup ##
#####################
## Drive field and sampling parameter

df_Amplitude = 20e-3    # 20 mT
df_Frequency = 26042
tx_Samples = 300        # 300 samples per period
magneticField =df_Amplitude*cos.(2*pi*(collect(1:tx_Samples)/tx_Samples))
timePoints = 1/df_Frequency * collect(1:tx_Samples)/tx_Samples

## Calibration coil parameter

R_cal = 100             # 100 Ohm
R_i = 50                # 50 Ohm
R = R_i*R_cal/(R_i+R_cal)
N = 15
r_cal = 1.25e-3
a = pi*r_cal^2
inductionFactor = R / (N*a)

#####################
## Generate MPS #####
# transfer function #

# Load transfer function
@info "Load S parameter"

S21_in_cal_a = TransferFunction(joinpath(folderName,"S21_in_cal_a.s1p"))
S21_in_cal_b = TransferFunction(joinpath(folderName,"S21_in_cal_b.s1p"))
S21_in_r = TransferFunction(joinpath(folderName,"S21_in_r.s1p"))

# Combine transfer function
S21_in_cal = vcat(S21_in_cal_a.data, S21_in_cal_b.data[length(S21_in_cal_a.data)+1:end]);
S21_in_cal = TransferFunction(S21_in_r.freq,S21_in_cal,[inductionFactor])

# Save MPS transfer function from u_meas to dm/dt
#
# We suggest this domain for MPS evaluation, since
# relaxation behavior can directly be seen.
# If the magnetic moment domain is desired one can
# simply ommit the derivation.
@info "Generate MPS transfer function"

derivationFactor = (2*pi.*S21_in_cal.freq*im)
G_dmdt = S21_in_cal.data./S21_in_r.data ./ derivationFactor;
G_dmdt = TransferFunction(S21_in_cal.freq,G_dmdt,[inductionFactor])

# Sample MPS transfer function to correct a shift of 5 samples of the cosine drive-field excitation
filename = joinpath(tmpdir,"example.mdf")
f = MPIFile(filename)
df_shift = -5
rxFrequs = rxFrequencies(f)
freqIndex = collect(1:length(rxFrequs))
G_dmdt_sampled = sampleTF(G_dmdt,f)
#close(f)
G_dmdt_sampled = G_dmdt_sampled.*exp.(df_shift/tx_Samples*2*pi*im.*collect(0:length(G_dmdt_sampled)-1))
G_dmdt_sampled = TransferFunction(rxFrequs,G_dmdt_sampled,[inductionFactor])

# if TF file already exists, delete it
if isfile("MPSTransferFunction_dmdt.h5")
  rm("MPSTransferFunction_dmdt.h5")
end
MPIFiles.save("MPSTransferFunction_dmdt.h5",G_dmdt_sampled)

@info "Set MPS transfer function"
setTF(f,"MPSTransferFunction_dmdt.h5")


#####################
### Visualization ###
#####################

# Load data in Frequency domain
@info "Load data"

u_meas_FD = squeeze(mean(getMeasurementsFD(MPIFile(filename),bgCorrection=true,tfCorrection=false),dims=4));
u_meas_FD = u_meas_FD./length(u_meas_FD)
dmdt_FD = squeeze(mean(getMeasurementsFD(MPIFile(filename),bgCorrection=true,tfCorrection=true),dims=4));
dmdt_FD = dmdt_FD./length(dmdt_FD)
m_FD = vcat(0,dmdt_FD[2:end]./(2*pi.*rxFrequs[2:end]*im)) 

# Load data in Time domain
u_meas_TD = squeeze(mean(getMeasurements(MPIFile(filename),bgCorrection=true,tfCorrection=false),dims=4));
dmdt_TD = squeeze(mean(getMeasurements(MPIFile(filename),bgCorrection=true,tfCorrection=true),dims=4));
m_TD = irfft(m_FD,300).*length(m_FD);

oddHarmonics = collect(2:2:151)
oddFreqs = rxFrequs[oddHarmonics]
evenHarmonics = collect(3:2:151)
evenFreqs = rxFrequs[evenHarmonics]

# absolute value of odd harmonics
u_meas_abs_odd = abs.(u_meas_FD[oddHarmonics])     
dmdt_abs_odd = abs.(dmdt_FD[oddHarmonics])  
m_abs_odd = abs.(m_FD[oddHarmonics]) 
# absolute value of even harmonics
u_meas_abs_even = abs.(u_meas_FD[evenHarmonics]) 
dmdt_abs_even = abs.(dmdt_FD[evenHarmonics]) 
m_abs_even = abs.(m_FD[evenHarmonics]) 

# angle of odd harmonics with corrections for 1exp(0deg)=-1exp(180deg) 
u_meas_angle_odd = angle.(u_meas_FD[oddHarmonics]).*180/pi     # in degree      
dmdt_angle_odd = angle.(dmdt_FD[oddHarmonics]) .*180/pi
dmdt_angle_odd[2:2:end] = angle.(-dmdt_FD[oddHarmonics[2:2:end]]) .*180/pi
m_angle_odd = angle.(m_FD[oddHarmonics]).*180/pi
m_angle_odd[2:2:end] = angle.(-m_FD[oddHarmonics[2:2:end]]).*180/pi

# angle of odd harmonics with corrections for 1exp(0deg)=-1exp(180deg) 
u_meas_angle_even = angle.(u_meas_FD[evenHarmonics]).*180/pi    
dmdt_angle_even = angle.(dmdt_FD[evenHarmonics]).*180/pi     
dmdt_angle_even[2:2:end] = angle.(-dmdt_FD[evenHarmonics[2:2:end]]).*180/pi
m_angle_even = angle.(-m_FD[evenHarmonics]).*180/pi      
m_angle_even[2:2:end] = angle.(m_FD[evenHarmonics[2:2:end]]).*180/pi  ;

# Create plots

# Time domain
p1 = plot(timePoints,u_meas_TD, xlabel = "t / µs", xticks = ([0,1e-5,2e-5,3e-5,4e-5],[0,1,2,3,4]), ylabel = "u_meas / V",lw = 1.3,legend = false)
p2 = plot(timePoints,dmdt_TD, xlabel = "t / µs", xticks = ([0,1e-5,2e-5,3e-5,4e-5],[0,1,2,3,4]), ylabel = "dm/dt / Am²/s",lw = 1.3,legend = false,title="Time domain")
p3 = plot(timePoints,m_TD, xlabel = "t / µs", ylabel = "m / µAm²", xticks = ([0,1e-5,2e-5,3e-5,4e-5],[0,1,2,3,4]), yticks = ([-6e-6,-4e-6,-2e-6,-0,2e-6,4e-6,6e-6],[-6,-4,-2,0,2,4,6]),lw = 1.3,legend = false)

# Frequency domain

# u_meas
# absolute
p4 = plot(oddFreqs,u_meas_abs_odd, yaxis=:log, xlabel = "f / MHz", ylabel = "|u_meas| / V",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([1e-8,1e-6,1e-4,1e-2,1],[1e-8,1e-6,1e-4,1e-2,1]),lw = 2,c="black",legend = false)
p4 = plot!(evenFreqs,u_meas_abs_even, yaxis=:log,lw = 2, c="darkorange2",legend = false)
# phase
p4 = plot!(twinx(),oddFreqs,u_meas_angle_odd,line=(:dash, 1.3),  c="black",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([-180,-135,-90,-45,0,45,90,135,180],[-180,-135,-90,-45,0,45,90,135,180]),ylim=(-180,180),legend = false)
p4 = plot!(twinx(),evenFreqs,u_meas_angle_even,line=(:dash, 1.3), c="darkorange2",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([-180,-135,-90,-45,0,45,90,135,180],[-180,-135,-90,-45,0,45,90,135,180]),ylim=(-180,180),legend = false,ylabel = "φ(u_meas) / °")

# m'
# absolute
p5 = plot(oddFreqs,dmdt_abs_odd, yaxis=:log, xlabel = "f / MHz", ylabel = "|m'| / Am²/s",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([1e-8,1e-6,1e-4,1e-2,1],[1e-8,1e-6,1e-4,1e-2,1]),lw = 2,c="black",label ="|signal_odd|")
p5 = plot!(evenFreqs,dmdt_abs_even, yaxis=:log,lw = 2,c="darkorange2",label ="|signal_even|",title="Frequency domain")
# phase
p5 = plot!(twinx(),oddFreqs,dmdt_angle_odd,line=(:dash, 1.3),  c="black",yticks = ([-180,-135,-90,-45,0,45,90,135,180],[-180,-135,-90,-45,0,45,90,135,180]),ylim=(-180,180),legend = false)
p5 = plot!(twinx(),evenFreqs,dmdt_angle_even,line=(:dash, 1.3),  c="darkorange2",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([-180,-135,-90,-45,0,45,90,135,180],[-180,-135,-90,-45,0,45,90,135,180]),ylim=(-180,180),legend = false,ylabel = "φ(m') / °")

# m 
# absolute
p6 = plot(oddFreqs,m_abs_odd, yaxis=:log, xlabel = "f / MHz", ylabel = "|m| / Am²",legend = false,xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([1e-13,1e-11,1e-9,1e-7,1e-5],[1e-13,1e-11,1e-9,1e-7,1e-5]),lw = 2,c="black")
p6 = plot!(evenFreqs[2:end],m_abs_even[2:end], yaxis=:log,lw = 2,c="darkorange2",legend = false)
# phase
p6 = plot!(twinx(),oddFreqs,m_angle_odd,line=(:dash, 1.3),  c="black",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([-180,-135,-90,-45,0,45,90,135,180],[-180,-135,-90,-45,0,45,90,135,180]),ylim=(-180,180),label ="φ(signal_odd)", fg_legend = :transparent)
p6 = plot!(twinx(),evenFreqs,m_angle_even,line=(:dash, 1.3),  c="darkorange2",xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),yticks = ([-180,-135,-90,-45,0,45,90,135,180],[-180,-135,-90,-45,0,45,90,135,180]),ylim=(-180,180),label ="φ(signal_even)", legend=:bottomleft, fg_legend = :transparent,ylabel = "φ(m) / °")

#Point Spread function
p7 = plot(magneticField,abs.(dmdt_TD),xticks = ([-0.02,-0.01,0,0.01,0.02],[-20,-10,0,10,20]),xlabel = "H / mT/µ₀", ylabel = "dm/dH / Am²/s",lw=1.3,title="Point Spread Function",legend = false)
#Hysteresis Curve
p8 = plot(magneticField,m_TD,xticks = ([-0.02,-0.01,0,0.01,0.02],[-20,-10,0,10,20]),xlabel = "H / mT/µ₀", ylabel = "m / Am²",lw=1.3,title="Hysteresis Curve",legend = false)

plot(p1,p2,p3,p4,p5,p6,p7,p8, layout = (3, 3), size=(1000,600),   
     #xticks = ([0,1e6,2e6,3e6,4e6],[0,1,2,3,4]),
     left_margin = 5mm,
     tickfontsize = 10)

