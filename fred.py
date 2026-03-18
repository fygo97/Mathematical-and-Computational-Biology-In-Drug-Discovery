'''
About the paper I found the following points interesting:

1. Vast Unexplored Space (this one I actually found absolutely mind blowing):  A significant portion of the genetically supported drug discovery space remains clinically unexplored. The researchers estimated that only 1.1% of all genetically supported gene-indication relationships have been explored clinically. Why? How do you make an educated decision that Gen A makes sense to explore and others not so much? How does this confidence align with the fact that roughly only 7% if all drugs make it through developments. Me as a mathematician would wonder "doesn't that mean that we are just completely off if our models yield that low a percentage and we should go back to fundamental research?". I do understand that 85% of the human body are considered "undruggable" but that still leaves about 13.9% of genes open to exploration as potential targets.

2. Pleiotropy: As a further explanation for the vast unexplored spaces I found that many of the 13.9% genes are hard to touch without starting a cascade of side effects. Why don't we focus more on fundamental research that allows us to learn how to target these genes safely but instead exhaust the drug targets that are already known and used?

3. Oncology paradoxon: Oncology is among those pathologies that has a large genetic component. How are drugs with a genetic support only negligibly more likely to succeed than drugs without genetic support? (compare figure 2.2) It would have made sense if oncologic drugs are just more successful overall but looking at the relative numbers of success this does not seem to be the case.


'''


import numpy as np
import matplotlib.pyplot as plt


# specificity = P(negative|healthy)  = 99.5% = 0.995
# sensitivity = P(positive|infected) = 99.5% = 0.995
# prevalence  = 5% = 0.05
#
# Since Fred took a test that came back positive we have 
# P(infected|positive) = P(infected and positive)/P(positive) by Bayes
#
# P(positive) = P(healthy and positive) + P(infected and positive)
#             = P(healthy) * P(positive|healthy)  + P(infected) * P(positive|infected)
#             = (1-prevalence)*(1-specificity) + prevalence*sensitivity
# => P(positive) = (1-0.05)*(1-0.995) + 0.05*0.995 = 0.00475+0.04975=0.0545 =:A
#  
# Thus 
#       P(infected|positive) = P(infected and positive)/P(positive) 
#                            = P(infected) * P(positive|infected)/ A
#                            = prevalence * sensitivity/ A
#                            = 0.995*0.05/0.0545
#                            = 91.28%
#
# which is "surprisingly low" but also shows why it is so very important to understand sensitivity and specificity
# when deploying tests as a medical professional (especially when handling very sensitive and rare diagnosis with a large implications; -> have seen 
# quite a few issues with this when I worked in the hospital as a nurse). Because that 91.28% means for a very common disease (5% prevalence is very high)
# still almost every 10th patient who might come in with:
#       - strong suspicion to have a certain disease
#       - fitting symptoms
#       - a positive test
# will be a false positive. Now extend this to something with a tiny prevalence like HIV with 0.19%
# (https://www.shcs.ch/news/kohler-et-al-hiv-care-cascade-in-switzerland/) and the risk of a false positive will increase from 10% to something much
# more common than 100-91.28%.


# Assumptions:
# positive test was taken
# specificity changes between 99%, 99,9%; 99,99%; 99,999%
# sensitivity fixed at 99%
# prevalence ranges from 0.001% to 50%




specificity = [0.99,0.999, 0.9999, 0.99999]
prevalence = np.linspace(0.00001, 0.5, 100000) # prevalence vector
P_positive = np.zeros((len(specificity),len(prevalence))) # the 
P_posterior = np.zeros((len(specificity),len(prevalence))) # the posterior P(infected|positive)

for i in range(len(specificity)):
    for j in range(len(prevalence)):

        #P(positive) = (1-prevalence)*(1-specificity) + prevalence*sensitivity
        P_positive[i,j] = (1-prevalence[j])*(1-specificity[i]) + prevalence[j]*0.99

        #P(infected|positive) = P(infected) * P(positive|infected)/ P(positive) = prevalence * sensitivity/ P(positive)
        P_posterior[i,j] = prevalence[j] * 0.99/ P_positive[i,j]

# Task 2: Visualization
plt.figure(figsize=(10, 6))
for i in range(len(specificity)):
    plt.plot(prevalence * 100, P_posterior[i] * 100, 
             label=f'Specificity: {round(specificity[i]*100, 3)}%')

plt.xscale('log') # Log scale makes the 0.001% range visible
plt.xlabel('Infection Prevalence (%)')
plt.ylabel('Prob. of being Infected if Positive (%)')
plt.title('PPV vs. Prevalence for Different Specificities')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.savefig('infection_probability_plot.png')
plt.show()

# Exercise 3
# The graph illustrates what I already teasered in line 31 of this code: Even for "nearly perfect" tests w.r.t. specificity (like 99.999%)
# the probability for a false positive essentially remains at the level of a coin flip for tiny prevalence (compare prevalence 0.00001,
# specificity 99.999% gives P(infected|positive) ~= 0.5 in the graph).
# It is also important to notice, that the graphs of the PPV (Positive Predictive Value) formula describe a sigmoid or hyperbolic shape, which means:
# For different levels of prevalence tests need to pass a certain threshold before a positive test becomes medically meaningful: For a test with 
# specificity of 99.999% (so an almost perfect test) crosses the P(infected|positive) = 90% threshold for prevalence of 10^-2 %=0.01% while a test with
# specificity 99.99% crosses the same probability threshold of P(infected|positive) = 90% for prevalence 0.1% (which is of magnitude 10 higher). So the 
# loss of safety in tests against false positive is fast as specificity deteriorates just a negligible amount.


