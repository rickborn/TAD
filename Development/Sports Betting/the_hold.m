function [h] = the_hold(odds,units)

% the_hold.m: calculate the sports book's margin
%
% h = the_hold(odds,units);
%
% e.g. h = the_hold([+130,-150],'American');
%
% Inputs:
% - odds, the odds offered for the two teams
% - units, what type of odds: American, Decimal or Fractional
%
% Outputs:
% - h, the hold percentage. Positive hold favors the book maker; negative
% hold favors the bettor.
%
% RTB wrote it, 07 February 2023, listening to Shostakovich string quartets
% Based on Miller & Davidow, "The Logic of Sports Betting" (2019)

% Calculate the break_even percentage for each of the odds
be = break_even(odds,units);

% Example: Imagine there's a Giants-Dodgers game, and the listed prices
% have the odds: Giants +130 and Dodgers -150. This means that the Giants
% are the underdog, and a $100 bet on them, if they win, wins you $130.
% Conversely, the Dodgers are the favorite, so you'd need to bet $150 on
% them in order to win $100 (if they win). The corresponding break-even
% percentages are 43.5% (Giants) and 60% (Dodgers). [This means that, over a
% large number of such bets, the Giants would need to win 43.5% of the time
% in order for you to break even. Similarly, the Dodgers would need to win
% 60% of the time--but remember that's why their the favorites.]
% 
% Now consider the situation in which you bet the Giants. If you buy the
% Giants at 43.5%, the sports book is "buying" the Dodgers at 56.5%. But note
% that the sports book is "selling" the Dodgers at 60%. It's this gap between
% what the sports book will buy a bet for and what he will sell it for that
% is known as "the hold." In this case, it's 60% - 56.5% = 3.5%. This is
% essentially the built-in profit margin of the book maker. Of course, this
% assumes equal betting on both sides, which the book maker tries to ensure
% by "price discovery," a process in which they post an opening line with
% relatively low betting limits, and then move the line according to how
% much action they take on either side. If all of the initial action picks
% the Giants, then they move the line to make them less of an underdog.
% 
% NOTE that we could make the calculation in the opposite direction. If you
% bet the Dodgers at 60%, it means the sports book is willing to buy the
% Giants at 40%, but only willing to sell them to you at a worse break-even
% percentage of 43.5%. Again, the difference (= "the hold") is 3.5%.
% 
% from pp. 23-4 of Miller & Davidow, "The Logic of Sports Betting" (2019)

% The math is pretty simple. It gets more complicated if one allows for
% ties (called a "push" in sports betting lingo), in which case everyone
% gets his money back, as if nothing ever happened.
h = be(2) - (100 - be(1));
