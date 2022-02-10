using(StatsBase)
using(Distributions)

function gillespie_wasps_fast(indata, t_max, param)
    t=0; i=1; j=1; alpha = param[1]; gamma = param[2]; beta = param[3]; tolerance=param[4]; time_switch_off=0;
    wasps = indata
    l=length(wasps)
    rates = zeros(1, 5*l)
    interactions = rand(Exponential(gamma+time_switch_off),1,10)
    ovaries = rand(Normal(0.1,0.005),1,10)
    gene_on = zeros(1, 10)
    pdec = zeros(1,l)
    counter = -100*ones(10)
    while t<t_max
    avg = sum(wasps).-wasps;
    pdec = cdf.(Normal(gamma+time_switch_off, 1e-5), interactions)
    rates = [500.*gene_on wasps alpha.*avg.*wasps beta.*pdec.*(1-gene_on)];
    dt = rand(Exponential(1/sum(rates)));
    t = t + dt;
    cumulative_expression = cumsum(rates,2)/sum(rates);
    rnumb = rand()
    rk = findfirst(x -> (x - rnumb) > 0, cumulative_expression)
    #Production
    if floor((rk - 1) / l) == 0
        wasps[mod(rk - 1, l) + 1] = wasps[mod(rk - 1, l) + 1] + 1.0;
    #Degradation
    elseif floor((rk - 1) / l) == 1
        wasps[mod(rk - 1, l) + 1] = wasps[mod(rk - 1, l) + 1] - 1.0;
    #Interactions - Gene flipping from 1 to 0
    elseif floor((rk - 1) / l) == 2
        attacked_wasp = mod(rk-1,l)+1
        attack_prob = copy(wasps)
        attack_prob[attacked_wasp] = 0;
        attack_rate = cumsum(attack_prob,2)/sum(attack_prob);
        atk_rand = rand();
        attacker = findfirst(x -> (x - atk_rand) > 0, attack_rate)
        if tolerance*(wasps[attacker]-wasps[attacked_wasp]) > 600
            sigmoidal_value = 1.0.*(wasps[attacker].>wasps[attacked_wasp])
        else
            sigmoidal_value = exp(tolerance*(wasps[attacker]-wasps[attacked_wasp]))/(exp(tolerance*(wasps[attacker]-wasps[attacked_wasp]))+1)
        end
        prob = rand()
        if rand() < sigmoidal_value # && interactions[attacked_wasp] > time_switch_off && gene_on[attacked_wasp] == 1
            interactions[attacked_wasp] = 0;
            gene_on[attacked_wasp] = 0;
        end
    #Gene (q) flipping from 0 to 1
    elseif floor((rk - 1) / l) == 3
        selected_wasp = mod(rk-1,l)+1
        gene_on[selected_wasp] = 1
    end
    wasps[find(x -> x<0, wasps)] .= 0;
    interactions = interactions .+ dt;
    gene_on[find(x -> abs.(x+time_switch_off-t)<dt, counter)].=0
    counter[find(x -> abs.(x+time_switch_off-t)<dt, counter)].=-100
    i = i+1
    end
    return wasps
end


sample_run = gillespie_wasps_fast(rand(Normal(50,15),1,10), 40, [1 1 1e5 10]);
