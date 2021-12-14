module fug

using Plots
using JuMP
using Ipopt

#geometric mean of relative volitilies, useful for shortcut column
"""
geoMean(kLKd, kLKb, kHKd, kHKb)

returns the geometric mean of k values

'''math
\\sqrt{LK_d*LK_b}{HK_d*HK_b}
'''

#Arguments 
*'kLKd': k value of light key in distilate
*'kLKb': k value of light key in bottoms
*'kHKd': k value of heavy key in distilate
*'kHKb': k value of heavy key in bottoms
"""
function geoMean(kLKd, kLKb, kHKd, kHKb)
    return(sqrt((kLKd*kLKb)/(kHKd*kHKb)))
end

#Fenske equation used in step 3 of column sizing to get the
"""
fenske(xLKd, kLKd, xHKd, kHKd, xLKb, kLKb, xHKb, kHKb)

returns the theoretical minimun number of trays for the distillation column

'''math
n = \\frac{log_{10}(\\frac{xLKd cdot xHKb}{xLKb cdot xHKd})}{log_{10}(\\alpha_{avg})}
'''


#Arguments 
*'xLKd': x value of light key in distilate
*'kLKd': k value of light key in distilate
*'xHKd': x value of heavy key in distilate
*'kHKd': k value of heavy key in distilate
*'xLKb': x value of light key in bottoms
*'kLKb': k value of light key in bottoms
*'xHKb': x value of heavy key in bottoms
*'kHKb': k value of heavy key in bottoms
"""
function fenske(xLKd, kLKd, xHKd, kHKd, xLKb, kLKb, xHKb, kHKb)
    α_avg = geoMean(kLKd, kLKb, kHKd, kHKb)
    n = log10(xLKd*xHKb/(xLKb*xHKd)) / log10(α_avg)
    return(n) #return number of trays
end

#nonKeyComp finds the compistions of the non-key components in both the distillate and bottoms
"""
nonKeyComp(alpha, nmin, F, D, B, xHKd, xHKb, xfi)

returns the molar fraction of a non-key component in both the bottoms and the distillate

(xBi, xDi)

by solving the 2 equations at the same time:
'''math
\\frac{x_{di}}{x_{bi}} = \\alpha^{n_{min}} \\cdot \\frac{x_{HK,d}}{x_HK,b}
'''

'''math
F_i = D_i + B_i
'''

done by using the following optimization:
'''math
Min[(\\frac{\\alpha^{N_{min}} \\cdot x_{HK,d}}{x_{HK,b}} - \\frac{x_{i,d} \\cdot B}{F \\cdot x_{i,f} - D \\cdot x_{i,d}})^2]
'''

#Arguments
*'alpha': relative volitility of component with the heavy key
*'nmin': the theoretical minimum number of stages for the column, get from Fenske equation
*'F': feed molar flow rate
*'D': distillate molar flow rate
*'B': bottoms molar flow rate
*'xHKd': x value of the heavy key in the distillate 
*'xHKb': x value of the heavy key in the bottoms
*'xfi': the x value of the non-key component in the feed
"""
function nonKeyComp(α, nmin, F, D, B, xHKd, xHKb, xfi)
    model = Model(with_optimizer(Ipopt.Optimizer))
    set_silent(model)
    @variable(model, xdi)
    @NLconstraint(model, xdi >= 1e-9)
    @NLconstraint(model, xdi <= xfi*F/D)
    @NLconstraint(model, (xfi*F-D*xdi)/B <= 1)
    @NLobjective(model, Min, (α^nmin*xHKd/xHKb - xdi*B/(F*xfi-D*xdi))^2)
    optimize!(model)
    xDi = value(xdi)
    xBi = (F*xfi-D*value(xdi))/B
    empty!(model)
    return(xBi, xDi)
end

#Underwood equation for finding theta
"""
underwood(kDist, kBot, kHKd, kHKb, xFeed, qual, xDist)

The Underwood equation for finding the minimum reflux ratio by solving the first equation subject to the second condition:

'''math
\\sum_{i=1}^n \\frac{\\alpha_i \\cdot x_{i,f}}{\\alpha_i - \\Theta} = 1-q
'''
'''math
R_{min} = \\sum_{i=1}^n \\frac{\\alpha_i \\cdot x_{i,d}}{\\alpha_i - \\Theta} - 1
'''

#Arguments
*'kDist': vector of k values for all components in the distillate
*'kBot': vector of k values for all components in the bottoms
*'kHKd': k value of heavy key in distillate
*'kHKb': k value of heavy key in bottoms
*'xFeed': vector of x values for the feed
*'qual': quality of the feed stream (percentage of vapor if in VLE)
*'xDist': vector of x values for the distillate
"""
function underwood(kDist, kBot, kHKd, kHKb, xFeed, qual, xDist)
    l = length(kDist)
    α = zeros(l)
    for i = 1:l
       α[i] = geoMean(kDist[i], kBot[i], kHKd, kHKb)
    end

    model = Model(with_optimizer(Ipopt.Optimizer))
    set_silent(model)
    @variable(model, Θ)
    @NLconstraint(model, 1.01 <= Θ <= maximum(α))
    @NLobjective(model, Min, (sum((α[i]*xFeed[i])/(α[i]-Θ) for i in 1:l) + qual -1)^2)
    optimize!(model)

    rmin = sum(((α .* xDist) ./ (α .- value(Θ)))) - 1
    empty!(model)
    return(rmin)
end

#Gilliand correlation for FUG method
"""
gilliand(ract, rmin, nmin)

The Gilliand correlation for the actual number of stages. Solve the following equation for N

'''math
\\frac{N-N_{min}}{N+1} = 0.75 \\cdot [1-(\\frac{R-R_{min}}{R+1})^{0.566}]
'''

#Arguments
*'ract': The actual reflux ratio, typically 1.1-1.4 times the minimum reflux ratio
*'rmin': Minimum reflux ratio (Underwood equation)
*'nmin': Theoretical minimum number of stages (Fenske equation)
"""
function gilliand(ract, rmin, nmin)
    model = Model(with_optimizer(Ipopt.Optimizer))
    set_silent(model)
    @variable(model, Nact)
    @NLconstraint(model, nmin <= Nact)
    @NLobjective(model, Min, (0.75*(1-((ract-rmin)/(ract+1))^0.566) - (Nact-nmin)/(Nact+1))^2)
    optimize!(model)
    nact = value(Nact)
    empty!(model)
    return(nact)
end

#Kirkbride function for finding optimal feed stage, uses ipopt solver
"""
kirkbride(B,D,xHKf,xLKf,xHKd,xLKb,nTot)

The Kirkbride function for finding the optimial feed stage. Solves the first equation subject to the constraint:

'''math
ln(\\frac{N_d}{N_b}) = 0.206 \\cdot ln(\\frac{B \\cdot x_{HK,f} \\cdot x_{LK,b}^2}{D \\cdot x_{LK,f} \\cdot x_{HK,d}^2})
'''
'''math
N_{tot} = N_d + N_b
'''

#Arguments
*'B': bottoms molar flow rate
*'D': distillate molar flow rate
*'xHKf': x value of heavy key in feed
*'xLKf': x value of light key in feed
*'xHKd': x value of heavy key in distillate
*'xLKb': x value of light key in bottoms
*'nTot': total number of stages in column
"""
function kirkbride(B,D,xHKf,xLKf,xHKd,xLKb,nTot)
    t = ((B/D) * (xHKf/xLKf) * (xLKb/xHKd)^2)^0.206
    model = Model(with_optimizer(Ipopt.Optimizer))
    set_silent(model)
    @variable(model, m) # where m is number of stages above feed
    @NLconstraint(model, 0.2 <= m <= nTot-0.2)
    @NLobjective(model, Min, (m/(nTot-m) - t)^2)
    optimize!(model)
    M = value(m)
    empty!(model)
    return(M, nTot-M) #return number stages above, number of stages below
end

"""
jj
"""
function diameter(F_ha, F_f, F_st, Cs, ρ_l, ρ_g, f, A_d, A_t, G)
    C = F_ha * F_f * F_st * Cs
    Uf = C*sqrt((ρ_l-ρ_g)/ρ_g)
    Dt = sqrt(4*G/ (f*Uf*pi*(1-A_d/A_t)*ρ_g))
    return(Dt)
end

### This is the master function for the fug method. The inputs will be:
### comp: List of component names
### xFeed: the component fractions in order from the least to the most volitile
### feedBasis: The molar flow rate of the feed
### kDist: the K values for the distillate
### kBot: the K values for the bottoms
### LKi: the index of the light key
### HKi: in the index of the heavy key
### splitFracLK: The percentage of the light key in the distillate for a 90/10 split this would be 0.9
### splitFracHK: the percentage of the heavy key in the bottoms
### F is molar flow rate of feed
### B is molar flow rate of bottoms
### D is molar flow rate of distillate
### quality is vapor fraction of the feed
### RRHeuristic is what to multiply the minimum reflux ratio by to get the actual reflux ratio


function fugk(comp, xFeed, kBot, kDist, HKi, LKi, splitFracHK, splitFracLK, F, B, D, quality,
    RRHeuristic)
    len = length(xFeed)
    feedMolFlow = xFeed * F
    botMolFlow=zeros(len); distMolFlow=zeros(len); xBot=zeros(len); xDist=zeros(len)
    
    botMolFlow[1:HKi-1]=feedMolFlow[1:HKi-1] 
    botMolFlow[HKi]=feedMolFlow[HKi]*splitFracHK
    botMolFlow[LKi]=feedMolFlow[LKi]*(1-splitFracLK)
    
    distMolFlow[HKi]=feedMolFlow[HKi]*(1-splitFracHK)
    distMolFlow[LKi]=feedMolFlow[LKi]*splitFracLK
    distMolFlow[LKi+1:end]=feedMolFlow[LKi+1:end]
        
    for i ∈ 1:len
        xDist[i] = distMolFlow[i]/sum(distMolFlow)
        xBot[i] = botMolFlow[i]/sum(botMolFlow)
    end
    print(xDist[HKi], " x val of HK in distilate " , xBot[LKi], " x val of LK in bottoms \n")
    
    xLKd = xDist[LKi]; kLKd = kDist[LKi]; xLKb = xBot[LKi]; kLKb = kBot[LKi];
    xHKd = xDist[HKi]; kHKd = kDist[HKi]; xHKb = xBot[HKi]; kHKb = kBot[HKi];

    nmin = fenske(xLKd, kLKd, xHKd, kHKd, xLKb, kLKb, xHKb, kHKb)
    print("Theoretical min number of stages is ", nmin, '\n')
    
    names=Array{Any}(undef,len-2); kNKd=zeros(len-2); kNKb=zeros(len-2); xfi=zeros(len-2);
    names[1:HKi-1] = comp[1:HKi-1]; names[HKi:end] = comp[LKi+1:end];
    kNKd[1:HKi-1] = kDist[1:HKi-1]; kNKd[HKi:end] = kDist[LKi+1:end];
    kNKb[1:HKi-1] = kBot[1:HKi-1]; kNKb[HKi:end] = kBot[LKi+1:end];
    xfi[1:HKi-1] = xFeed[1:HKi-1]; xfi[HKi:end] = xFeed[LKi+1:end];
    xHKd = xDist[HKi]
    xHKb = xBot[HKi]
    
    alp = zeros(len-2)
    for i = 1:len-2
         alp[i] = geoMean(kNKd[i], kNKb[i], kHKd, kHKb)
    end

    xDist=zeros(len-2); xBot=zeros(len-2);
    for i = 1:len-2
        xBot[i], xDist[i] = nonKeyComp(alp[i], nmin, F, D, B, xHKd, xHKb, xfi[i])
        print(names[i], " bottoms frac ", xBot[i], " distilate frac ", xDist[i], '\n')
    end
    
    insert!(xBot,HKi,xHKb); insert!(xBot,LKi,xLKb);
    insert!(xDist,HKi,xHKd); insert!(xDist,LKi,xLKd);
    
    rmin = underwood(kDist, kBot, kHKd, kHKb, xFeed, quality, xDist)
    print("The min reflux ratio is ", rmin, '\n')
    ract = rmin * RRHeuristic
    
    nact = gilliand(ract, rmin, nmin)
    print("Our actual number of stages will be ", nact, '\n')
    
    above, below = kirkbride(B,D,xFeed[HKi],xFeed[LKi],xLKb,xHKd, nact)
    print("our optimal feed stage is ", above, " from the top and ", below, " from the bottom \n")

end

end # module
