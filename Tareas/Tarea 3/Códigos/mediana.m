function [med] = mediana(EDP,IM)

asd = zeros(1,length(EDP));
new_IM = (0.1:0.1:20).';

for i = 1 : length(EDP)
    n = 0;
    mult = 1;
        for j = 1 : 20
            a = EDP(i,j);
            b = round(IM(i,j));
            c = round(new_IM(i));
            if ~isnan(a) 
                if  b == c
                    n = n + 1;
                    mult = mult * a;
                end
            end
        end
    asd(i) = mult^(1/n);
end
med = asd';
end

