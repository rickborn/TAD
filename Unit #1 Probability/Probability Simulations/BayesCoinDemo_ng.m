% BayesCoinDemo_ng.m
%
% Based on http://www.nature.com/nmeth/journal/v12/n5/full/nmeth.3368.html
%
% RTB wrote it 04 May 2017 (badly itching eyes)

% This exercise is designed to walk you through a simple probability
% problem from two different perspectives: frequentist vs. Bayesian. It is
% based on an excellent article in Nature Method's series on statistics,
% called 'Points of Significance: Bayesian statistics'
%
% Jorge L�pez Puga, Martin Krzywinski & Naomi Altman, Nat. Methods 12,
% 377�378 (2015);

% What to do: Login to learning catalytics and join the session for the
% modules entitled "Bayes coin demo" and "Bayes coin demo, Long answer."
% You will answer a series of questions based on the guided programming
% below. There are two types of questions. First, those that have a
% numerical or multiple choice answer will be done in the first module as a
% "Team-based assessment": The first time through, you will attempt to
% answer the questions on your own; the second time you will work with your
% team and submit one answer per team. Don't worry about your scores! (The
% class is pass/fail. This is just a way to have fun and track your
% progress.) The second type are those that are more philosophical and have
% longer answers. Use the 2nd module to answer these.
%
% Each section begins with a '%%'. Read through the comments and follow the
% instructions provided. In some cases you will be asked to answer a
% question, clearly indicated by 'QUESTION'. In other cases, you be asked
% to supply missing code, indicated by 'TODO'. The corresponding question
% in learning catalytics will be indicated in parentheses (e.g. Q1 for the
% 1st module; L1 for the 2nd). If there is no 'Q#' accompanying a
% 'QUESTION' just type your answer into this script and discuss it with
% your team. Once you have supplied the required code, you can execute that
% section by mouse-clicking in that section (The block will turn yellow.)
% and then simultaneously hitting the 'ctrl' and 'enter' keys (PC) or
% 'command' and 'enter' keys (Mac).

%% We flip a coin 3 times and get 3 heads: Frequentist approach
% We want to know whether the coin is fair or not.

% The frequentist approach is to ask how likely our data is under the null
% hypothesis (H0)
% QUESTION (Q1): What is H0 in terms of the probability of getting a head?

% TODO (Q2): Under H0, calculate the probability of getting 3 heads on 3 tosses
p3H = 

% QUESTION (Q3): What would our probability be if we allowed for either 3 heads
% OR 3 tails on 3 tosses?

% TODO (Q4): Use 'binofit' to calculate a 99% confidence interval for P(heads)
% given that we got 3 heads on 3 tosses.


% QUESTION (Q5): What happens to our confidence interval as we increase the
% number of tosses?

% QUESTION (Q6): What do we mean by a 'confidence interval':


%% We flip a coin 3 times and get 3 heads: Bayesian approach

% Here we consider an entire range of possible hypotheses concerning the
% probability of heads. That is, we consider the full range from 
% P(heads) = 0 to P(heads) = 1.
pi0 = 0:0.001:1;
nHypotheses = length(pi0);

% NOTE: It's important to fully grok what our 'pi0' variable is. It is a
% vector of 1001 values, each of which should be thought of as a hypothesis
% about the true value of p(Heads) for our coin. For example, our 501st
% hypothesis is that the coin is fair, i.e. pi0(501) = 0.5, but this is
% just one of many that we are willing to consider. We start with some
% prior set of beliefs about all of the possible values of p(Heads). What
% Bayesian inference will allow us to do is, based on observed data,
% reallocate credibility across the different possibilities.

% If we have no clue as to what the true value of P(heads) is, we could
% consider all the values to be equally likely. This is often referred to
% as a 'flat prior' (because the histogram is flat) or an 'uninformative
% prior' (because all of our hypotheses about pi0 are equally likely).
priorFlat = ones(size(pi0)) ./ nHypotheses;

% QUESTION (L1): How intuitively plausible is this prior? Imagine that you have
% visually inspected the coin before we did our experiment, and you saw
% that, indeed, one side of the coin was a head and the other was a tail.

% Now, similarly to  the frequentist, we want to know the probability of
% our data (3 heads on 3 tosses) given a hypothesis. This is called the
% 'likelihood'. But here, as Bayesians, we are going to calculate the
% likelihood for each possible value of P(heads). That is, for all of the
% elements in our vector 'pi0'

% TODO: Calculate the probability of 3 heads on 3 tosses for each element
% in pi0 and store this in a variable called 'likelihood3H'. HINT: You can
% use a 'for' loop, but you don't have to.
likelihood3H = ;
likelihood3H = likelihood3H ./ trapz(likelihood3H); % normalize by area

% NOTE: We don't need to normalize the likelihood to do a proper
% calculation of our posterior, since we are going to normalize the
% posterior anyway. We do it here for plotting purposes.

% Note that the frequentist p-value ('p3H' calculated above) is a subset of
% the likelihood we just calculated. That is: likelihood3H(pi0 == 0.5)

% Now we use Bayes' rule to calculate the probability of each hypothesis
% (i.e. the possible values of P(heads) in pi0), known as the 'posterior
% probability': P(hypothesis|data). This is just the product of our prior
% and our likelihood, followed by normalization.
% TODO: Calculate the posterior distribution.
posterior3H = ;
posterior3H = posterior3H ./ trapz(posterior3H);    % normalize

% Plot our results
figure, 
subplot(2,2,1)
plot(pi0,priorFlat,'k-',pi0,likelihood3H,'r-',pi0,posterior3H,'b--');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Bayes: flat prior');

%% Calculate some values of interest

% TODO (Q7): Calculate the mean value for the posterior distribution.
% HINT: Recall that the mean value is the expected value. You can't just
% take the mean of the posterior, since it is a probability distribution,
% not a sample from that distribution. This is a subtle but important
% concept. If it is not clear to you, ask one of the instructors.
pi0mean = ;

% TODO (Q8): Calculate the median value for the posterior distribution.
% HINT: The median divides the probability mass in half.
pi0median = ;

% TODO (Q9): Calculate the most likely value of pi0 (i.e. the mode):
% HINT: Despite the similar names, this is not about the likelihood!
% Think about what the posterior distribution is really saying: this is the
% probability of each value of pi0, given what we already think we know
% about pi0s in general.
pi0max = ;

% But, in addition, the Bayesian can answer interesting questions that 
% don't even exist for the frequentist.

% TODO (Q10): Calculate the probability that the coin has any amount of bias towards heads?
pBiased = ;

% QUESTION (Q11) What are the odds the coin is biased towards heads?
oddsBiased2Heads = pBiased / (1-pBiased);

% QUESTION (Q12): How do we interpret this odds value? If you were going to bet
% 1$ that the coin was fair (i.e. NOT biased towards heads), how many $
% should the pay-out be in order to make it a fair bet?
% HINT: see http://andrewgelman.com/2010/07/10/creating_a_good/

%% Non-flat priors
% Now, what if we base our prior on having previously observed 3H and 1T:

% TODO: Calculate a prior based on the observation of 3 heads and 1 tails
% HINT: If we originally thought the distribution was uniform (pi0), then
% another way to think about this prior is as the posterior distribution
% after we've already observed those 4 trials.
% (For plotting purposes, re-normalize each of your variables.)
prior3H1T = ;

% Then we observe an outcome of 4H:
likelihood4H = binopdf(4,4,pi0);

% TODO: Compute the posterior, using our likelihood and prior
posterior3H1T = ;

% Now plot the new posterior distribution
subplot(2,2,2);
plot(pi0,prior3H1T,'k-',pi0,likelihood4H,'r-',pi0,posterior3H1T,'b-');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Prior based on 3H1T');

% QUESTION (L2):
% What is the difference between a likelihood and the posterior?
% How is it that Bayes rule is changing what we think about pi0?

% QUESTION (L3):
% Is this result any different from using the uniform prior and then
% observing 7 heads and 1 tail? Should it be?

%% Gaussian prior centered around 0.5, observe 3H1T

% Now let's consider a slightly more realistic situation. Since we've had a
% chance to examine our coin and we see that it looks like a relatively
% normal disc with 'heads' on one side and 'tails' on the other, we are
% reasonably sure, a priori, that it is a fair coin. However, we are still
% willing to entertain the notion that the coin might be biased in one
% direction or the other. We can implement this prior using a normal
% distriubution centered around 'fairness' (i.e. p(H) = 0.5) but having a
% width that reflects our uncertainty about the fairness hypothesis:
priorGaussFair = normpdf(pi0,0.5,0.25);
priorGaussFair = priorGaussFair ./ trapz(priorGaussFair);

% Then we observe an outcome of 3H1T:
% TODO: Compute the likelihood and posterior distributions
likelihood3H1T = ;

posterior3H1T = ;

% QUESTION (Q13): Given this posterior, calculate the probability that the
% coin has any amount of bias towards heads?

subplot(2,2,3);
plot(pi0,priorGaussFair,'k-',pi0,likelihood3H1T,'r-',pi0,posterior3H1T,'b-');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Gaussian prior around 0.5');

% QUESTION (L4):
% How is this prior different from the one we imposed above using "previous
% observations?" Not just numerically, but philosophically. What is the ideological
% distinction between the two types of priors?
%% Gaussian prior tightly centered around 0.5, observe 3H1T

% Now we might imagine that we have some data from the US Government
% indicating that most coins are very close to being fair. That is, we are
% not very willing to entertain hypotheses like p(H) = 1.

% TODO: Create a Gaussian prior with a standard deviation of 0.1
priorGaussFairNew = ;

% Now we flip the coin 4 times and observe, 3H,1T:
likelihood3H1TNew = ;

posterior3H1TNew = ;

subplot(2,2,4);
plot(pi0,priorGaussFairNew,'k-',pi0,likelihood3H1TNew,'r-',pi0,posterior3H1TNew,'b-');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Prior tight near 0.5');

% Note that in both of our bottom subplots, the likelihood (red curve) is
% identical but our posterios appear to have shifted by different amounts.
% To what do you attribute this? Does it make sense?

% TODO: As a measure of the shifts, compare your assessment of the
% probability that the coin has any amount of bias towards heads:
pBiasedNowNew = ;

% TODO (L5): Save your final figure as a jpeg, then upload it to the
% learning catalytics module.