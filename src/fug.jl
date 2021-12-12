module fug

using Plots
using JuMP
using Ipopt

#geometric mean of relative volitilies, useful for shortcut column
function geoMean(kLightKeyDist, kLightKeyBot, kHeavyKeyDist, kHeavyKeyBot)
    return(sqrt((kLightKeyDist*kLightKeyBot)/(kHeavyKeyDist*kHeavyKeyBot)))
end

#Fenske equation used in step 3 of column sizing to get the
function fenske(xLightKeyDist, kLightKeyDist, xHeavyKeyDist, kHeavyKeyDist, xLightKeyBot, kLightKeyBot,
        xHeavyKeyBot, kHeavyKeyBot)
    α_avg = geoMean(kLightKeyDist, kLightKeyBot, kHeavyKeyDist, kHeavyKeyBot)
    n = log10(xLightKeyDist*xHeavyKeyBot/(xLightKeyBot*xHeavyKeyDist)) / log10(α_avg)
    return(n) #return number of trays
end

function step4(α, nmin, F, D, B, xHKd, xHKb, xfi)
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
