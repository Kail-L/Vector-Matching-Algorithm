function dyear = decYear_mod(year, month, day)
%DECYEAR_MOD  Converts day/month/year date into decimal year date.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DYEAR = decYear_mod(year, month, day) This function converts date given
% in day/month/year into a decimal year. It does not account for leap year
% values, which can be implemented in future iterations.
% 
% SOURCES:
% n/a (created by previous ADCS team member)
%
% INPUT PARAMETERS:
% year = integer year of choice (ex. 2019)
% month = integer month of choice (ex. 5 corresponding to May)
% day = integer day of choice (ex. 23)
%
% OUTPUT PARAMETERS:
% dyear = integer value of decimal year
%
% VARIABLES:
% monthLengths = array containing number of days in given month.
%              = [Jan. Feb. Mar...]
% monthDays = number of days up to given input date
% daysInYear = integer number of days in given year
%
% Kail Laughlin, taken from uknown previous ADCS team member
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

monthLengths = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
if month > 1
    monthDays = sum(monthLengths(1:month-1));
else
    monthDays = 0;
end
daysInYear = 365.25;

dyear = year + (day + monthDays)/daysInYear;
end