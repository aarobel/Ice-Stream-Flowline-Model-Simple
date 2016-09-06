# Ice-Stream-Flowline-Model-Simple

# This is a numerical implementation of the flowline model described in Schoof, JGR, 2007 for MATLAB. It solves for ice thickness, grounding line position and velocity, given a bed topography (defined in Base.m and dBasedx.m), an accumulation rate, Glen's Law Rate Factor and many other parameters. The default version of the model uses a Weertman-type power law sliding and bed topography from Schoof 2007. You can find a detailed mathematical description of the numerics in Robel et al., JGR, 2014. 

# Plug and play instructions:
# The model runs from the GroundingLine_FlowlineModel.m function. Inputs are accumulation rate and Glen's law rate factor and outputs are time-dependent grounding line positions, and icethickness for 10 kyr. In this sort of model, you can't easily initialize it from any arbitrary grounding line position, so I've created one initial condition (which is included) that corresponds to the large steady-state for accumulation rate of 0.3 m/yr and A_glen = 1e-25 . If you put in these parameters, it will stay in the same place, if you use other parameters, it will evolve away from that initial condition. Don't change A_glen too much though, because then it has trouble finding an initial solution.

# If you set the "plotting" input parameter to 1, it will also plot some things from the model as it runs. You can quickly test that it works by running it with something like: [time_all,xg_all,h_all,parameters] = GroundingLine_FlowlineModel(0.2/(3600*24*365),2e-25,1);

# I've set the model to be relatively coarse (~1 km resolution, 100 m near the grounding line) so that it runs fast on a laptop, and you can see the output in real time.
