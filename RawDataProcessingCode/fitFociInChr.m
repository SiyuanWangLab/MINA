function [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity, Xwidth, Ywidth, Zwidth] = fitFociInChr(ImageStack, IchrIndex, Zrange)
% updated on 190314.

% updated on 181203 to include median filtering before finding potential
% center postion, and then use the unfiltered image to fit the exact
% center. The reason for the change is that the previous strategy of
% fitting the center with filtered image may cause a minor shift in the
% fitted center position.

ImageMax = medfilt2(max(ImageStack, [],3));
S = size(ImageMax);
for ii = 1:length(IchrIndex)
    Ichr = zeros(size(ImageMax));
    Ichr(IchrIndex{ii}) = 1;
    ImageMaxChr = ImageMax.*Ichr;
    Ind = find(ImageMaxChr == max(max(ImageMaxChr)));
    [Y0, X0] = ind2sub(S, Ind(1));
    [C, I] = max(permute(ImageStack(Y0, X0,Zrange{ii}(1):Zrange{ii}(2)), [3 1 2]));% corrected on 190606
    I = I-1+Zrange{ii}(1);% corrected on 190606
    CenterIntensity = ImageStack(Y0, X0, I);

    Data1 = mean(ImageStack(Y0,X0-7:X0+7,I),1);
    Data2 = mean(ImageStack(Y0-7:Y0+7,X0,I),2);

    GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
    StartPoint1 = [CenterIntensity X0 1 0];
    StartPoint2 = [CenterIntensity Y0 1 0];

    [f1, Xgof(ii)] = fit((X0-7:X0+7)', Data1', GaussEqu, 'Start', StartPoint1);
    [f2, Ygof(ii)] = fit((Y0-7:Y0+7)', Data2, GaussEqu, 'Start', StartPoint2);

    Zprofile = permute(ImageStack(Y0, X0, :), [3 1 2]);
    if I-3<1
        Irange = 1:7;
    elseif I+3> length(Zprofile)
        Irange = length(Zprofile)-6:length(Zprofile);
    else Irange = (I-3:I+3);
    end

    Data3 = Zprofile(Irange);
    StartPoint3 = [CenterIntensity I 1 0];
    [f3, Zgof(ii)] = fit(Irange', Data3, GaussEqu, 'Start', StartPoint3);

    if f1.b>X0-7 && f1.b<X0+7 && f2.b>Y0-7 && f2.b<Y0+7 && f3.b>Irange(1) && f3.b<Irange(end) %&& f1.c<2 && f2.c<2
        Xfit(ii) = f1.b;
        Yfit(ii) = f2.b;
        Zfit(ii) = f3.b;

        Intensity(ii) = (f1.a*f2.a*f3.a)^(1/3);
        Xwidth(ii) = f1.c;
        Ywidth(ii) = f2.c;
        Zwidth(ii) = f3.c;
    else
        Xfit(ii) = nan;
        Yfit(ii) = nan;
        Zfit(ii) = nan;
        Intensity(ii) = nan;
        Xwidth(ii) = nan;
        Ywidth(ii) = nan;
        Zwidth(ii) = nan;
    end
end
