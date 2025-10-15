function [be] = break_even(odds,units)

% break_even.m: calculate the break-even percentage based on odds
% 
% A break-even percentage (alternatively called "implied odds") is the
% percentage of time a bet must win for you to neither win nor lose money,
% when making the bet over time. If someone offers you an even money bet
% that a coin flip will land heads, the break-even percentage is 50%. Not
% because those are the true odds of winning, but because that's how often
% you have to win to exactly break even on an even money bet.
% 
% Break-even percentage = Risk / (Risk + Win)
% 
% from Miller & Davidow, "The Logic of Sports Betting" (2019)
% 
% Inputs:
% - odds, the odds offered for the game
% - units, what type of odds: American, Decimal or Fractional
% 
% Outputs:
% - be, the break-even percentage
% 
% RTB wrote it, 06 February 2023, post-RYP tutoring of Eva in trig
% 
% NOTE on different types of odds:
% 
% American odds are represented by either a positive number or a negative
% number, and the number is always at least 100. The minus sign denotes
% that the bet is on a favorite. The plus sign denotes an underdog. When
% the bet is plus, the number after the plus sign means how much money you
% would win if you bet $100 on that market. So if you make a bet at +350
% odds, it means you will win $350 for every $100 you bet. If you bet $10,
% you win $35. If you bet $200, you win $700. And so on. When the bet is
% minus (i.e., a favorite) the number after the minus sign means how much
% money you would have to bet to win $100 on that market. So if it’s -110
% that means you have to bet $110 if you want to win $100. An even-money
% bet would be +100 or -100 or EV.
% 
% Decimal odds are quoted as a single number, such as 1.9 or 2.5 or 6. What
% that number means is it’s the amount that you’ll collect from a winning
% $1 bet—including the amount of your bet. So if you bet $1 at 2.5, you
% will receive $2.50 from a winning bet— your original $1 plus $1.50 in
% winnings. Therefore, 2 represents an even money bet. Numbers between 1
% and 2 represent favorites. And numbers over 2 represent underdogs.
% 
% Fractional odds. If you see something like 7/4 odds, it means if you bet
% 4 (the denominator) you stand to win 7 (the numerator). In American, 7/4
% is +175. In decimal it’s 2.75. You would represent -110 as 10/11 in
% fractional odds. Even-money would be 1/1 odds.
% 
% It looks like the relevant parts of the book are available online:
% https://issuu.com/sbc.global/docs/the_logic_of_sports_betting_v1_0_5/s/12292345

switch units
    case 'American'
        disp('Using American odds.');
        be = zeros(size(odds));
        L = odds > 0;
        be(L) = 100 ./ (100 + odds(L));
        be(~L) = -(odds(~L)) ./ (100 + -(odds(~L)));
        
%         be = zeros(size(odds));
%         for k = 1:length(odds)
%             if odds(k) > 0
%                 be(k) = 100 ./ (100 + odds(k));
%             else
%                 pos_odds = -odds(k);
%                 be(k) = pos_odds ./ (100 + pos_odds);
%             end
%         end
    case 'Decimal'
        disp('Using Decimal odds.');
        be = 1 ./ odds;
    case 'Fractional'
        disp('Using Fractional odds.');
        be = 1 ./ (odds + 1);
    otherwise
        be = NaN;
        error('Units for odds incorrectly specified.');
end

% convert to a percentage:
be = round(be * 100,2);

end

