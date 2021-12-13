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
function underwood(kDist, kBot, kHeavyKeyDist, kHeavyKeyBot, xFeed, qual, xDist)
    l = length(kDist)
    α = zeros(l)
    for i = 1:l
       α[i] = geoMean(kDist[i], kBot[i], kHeavyKeyDist, kHeavyKeyBot)
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
function kirkbride(B,D,xHKFeed,xLKFeed,xLKBot,xHKDist, nTot)
    t = ((B/D) * (xHKFeed/xLKFeed) * (xLKBot/xHKDist)^2)^0.206
    model = Model(with_optimizer(Ipopt.Optimizer))
    set_silent(model)
    @variable(model, m) # where m is number of stages above feed
    @NLconstraint(model, 0.2 <= m <= nTot-0.2)
    @NLobjective(model, Min, (m/(nTot-m) - t)^2)
    optimize!(model)
    M = value(m)
    print("test")
    empty!(model)
    return(M, nTot-M) #return number stages above, number of stages below
end

end # module
