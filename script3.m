% script3.m
% use MATLAB's own pricing methods

% T = 1; Sx = 0.4; Sy = 0.4; rho = 0; Rf = 0.04; K = 100;
Settle = '01-Jan-2012';
Maturity = strcat(['01-Jan-',int2str(2012+T)]);
Price1 = K;
Vol1 = Sx;
Price2 = K;
Vol2 = Sy;
Corr = rho;
OptSpec = 'call';
Strike = 0;
rates = Rf;
Compounding = -1;
Basis = 1;
RateSpec = intenvset('ValuationDate', Settle, 'StartDates', Settle, ...
'EndDates', Maturity, 'Rates', rates, ...
'Compounding', Compounding, 'Basis', Basis);
StockSpec1 = stockspec(Vol1, Price1);
StockSpec2 = stockspec(Vol2, Price2);

Price = spreadbyls(RateSpec, StockSpec1, StockSpec2, Settle, ...
Maturity, OptSpec, Strike, Corr, 'AmericanOpt', OptionType);
disp(strcat(['MATLAB LSMC: ',num2str(Price)]));

Price = ...
 spreadbyfd(RateSpec, StockSpec1, StockSpec2, Settle, ...
 Maturity, OptSpec, Strike, Corr, 'AmericanOpt', OptionType);
 disp(strcat(['MATLAB ADI: ',num2str(Price)]));