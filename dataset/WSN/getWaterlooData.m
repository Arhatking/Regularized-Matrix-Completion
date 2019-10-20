function X = getWaterlooData
    load('waterloo_04.txt');
    load('waterloo_05.txt');
    load('waterloo_06.txt');
    load('waterloo_07.txt');
    load('waterloo_08.txt');
    load('waterloo_09.txt');
    X=[waterloo_04(1:96*360)' waterloo_05(1:96*300)' waterloo_06(1:96*360)' waterloo_07(1:96*360)' waterloo_08(1:96*360)' waterloo_09(1:96*360)'];
    X = reshape(X,96,360*5+300);
    
    X_avr = mean(X);
    idx = find(X_avr>20);
    X = X(:,idx);
end

