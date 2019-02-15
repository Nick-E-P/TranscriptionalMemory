function [yTOTranded] = Randomise_yTOT(yTOT)
%RANDOMISE_YTOT Randomises pairings of cells to control for similarity in
%cell-cycle length between sisters

for w = 1:size(yTOT,2)
division_length(w) = find(yTOT(:,w),1,'last')*5/60;
end

firstL = division_length(1:2:end);
secondL = division_length(2:2:end);

firstI = (1:2:size(yTOT,2));
secondI = 2:2:size(yTOT,2);

sigma = std((firstL-secondL));

secondLcurr = secondL;
secondIcurr = secondI;

for k = 1:round(size(yTOT,2)/2)
    trial = firstL(k)+normrnd(0,sigma);
    [~, idx] = min(abs(trial - secondLcurr));
    while length(secondLcurr)~=1 && secondIcurr(idx) == k+1
        trial = firstL(k)+normrnd(0,sigma);
        [~, idx] = min(abs(trial - secondLcurr));
    end
    listS(k) = secondIcurr(idx);
    secondLcurr(idx) = [];
    secondIcurr(idx) = [];
    
end

yTOTranded(:,1:2:size(yTOT,2)) = [yTOT(:,firstI)];
yTOTranded(:,2:2:size(yTOT,2)) = [yTOT(:,listS)];
end

