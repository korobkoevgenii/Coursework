using DifferentialEquations, FindPeaks
using Plots
my=1000
C=1
gNa=40
gK=35
gLeak=0.3
ENa=55
EK=-77
ELeak=-65
Iapp1=1.3
Iapp2=1.2


function neyron(du,u,p,t)
    #V,m,n,h,V2,m2,n2,h2=u
    V=u[1]
    m=u[2]
    n=u[3]
    h=u[4]
    
    V2=u[5]
    m2=u[6]
    n2=u[7]
    h2=u[8]
    # first neyron
    du[1]=my*(1/C)*(gNa*(m^3)*h*(ENa-V)+gK*n*(EK-V)+gLeak*(ELeak-V)+Iapp1)
    du[2]=my*(((0.182*(V+35))/(1-exp((-(V+35))/9)))*(1-m)-(((-0.124)*(V+35))/(1-exp((V+35)/9)))*m)
    du[3]=my*(((0.02*(V-25))/(1-exp((-(V-25))/9)))*(1-n)-(((-0.002)*(V-25))/(1-exp((V-25)/9)))*n)
    du[4]=my*((0.25*exp((-(V+90))/12))*(1-h)-(0.25*((exp((V+62)/6))/(exp((V+90)/12))))*h)
    
    # second neyron
    du[5]=my*(1/C)*(gNa*(m2^3)*h2*(ENa-V2)+gK*n2*(EK-V2)+gLeak*(ELeak-V2)+Iapp2)
    du[6]=my*(((0.182*(V2+35))/(1-exp((-(V2+35))/9)))*(1-m2)-(((-0.124)*(V2+35))/(1-exp((V2+35)/9)))*m2)
    du[7]=my*(((0.02*(V2-25))/(1-exp((-(V2-25))/9)))*(1-n2)-(((-0.002)*(V2-25))/(1-exp((V2-25)/9)))*n2)
    du[8]=my*((0.25*exp((-(V2+90))/12))*(1-h2)-(0.25*((exp((V2+62)/6))/(exp((V2+90)/12))))*h2)

end

# возвращаем массив по решению
# solve само решение
# comp
function getArrFromSol(solve,comp)
    
    size1=size(solve.t)[1]
    out=zeros(size1)
    if(comp==0)
        out=solve.t
    else
        for i in 1:size1
            out[i]=solve.u[i][comp]
        end
    end
    return out
end


# Функция возвращает локацию во времени по результатам поиска пиков
function GetTrueLocation(locations, ratio)
    trueLocation=Array{Float64}(undef,0)
    push!(trueLocation, 0.0)
    for loc in locations
        newLoc=loc/ratio
        push!(trueLocation, newLoc)
    end
    return trueLocation
end


# Функция возвращает частоту по местоположению пиков
function GetFrequency(locations)
    
    frequency=Array{Float64}(undef,0)
    prevLoc=0.0
    prevPeriod=0.0
    
    for location in locations
        period=location-prevLoc
        freq= prevPeriod==0.0 ? 0.0 : (1.0/period)
        prevPeriod=period
        push!(frequency, freq)
        prevLoc=location
    end
    
    return frequency
end

#u0 = [-58.7085;0.0953;0.000913;0.3662;-58.7085;0.0953;0.000913;0.3662] #stable focus
u0 = [14.8409;0.9174;0.0140;0.0539;14.8409;0.9174;0.0140;0.0539] #limit cycle
#u0 = [1.0;0.0;0.0;0.0]
tspan = (0.0,1)
prob = ODEProblem(neyron,u0,tspan)
sol = solve(prob,RK4())

time=getArrFromSol(sol,0)
firstV=getArrFromSol(sol,1)
secondV=getArrFromSol(sol,5)

plot(time, [firstV, secondV], label="V(t)", yaxis="V", xaxis="t")
#plot(sol, vars=(0,5), label="V(t)", yaxis="V", xaxis="t")


firstPeaks = findpeaks(firstV)
secondPeaks = findpeaks(secondV)
p1=plot(firstPeaks, firstV)
p2=plot(secondPeaks, secondV)
plot(p1,p2,layout = (2,1))

locationsFirst=peaklocations(firstPeaks)
heightsFirst=peakheights(firstPeaks)
locationsSecond=peaklocations(secondPeaks)

trueLocation1=GetTrueLocation(locationsFirst,4000)
trueLocation2=GetTrueLocation(locationsSecond,4000)




frequency1=GetFrequency(trueLocation1)
frequency2=GetFrequency(trueLocation2)

plot(trueLocation1, frequency1, linetype=:steppre, label="Freq1(t)", yaxis="Frequency", xaxis="t")
plot!(trueLocation2, frequency2, linetype=:steppre,  label="Freq2(t)")


