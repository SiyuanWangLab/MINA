function yy = generateSECDED16(NumOnes)
%generate SECDED(16,11) codes with the given number of ones

ii=0;
for i = 1:2^11-1
    BinString = dec2bin(i);
    N = length(BinString);
    DataBits = BinString(N:-1:1);
    if N<11
        for j = 1:11-N
            DataBits = [DataBits '0'];
        end
    end
    p1 = mod(str2num(DataBits(1))+str2num(DataBits(2))+str2num(DataBits(4))+str2num(DataBits(5))+str2num(DataBits(7))+str2num(DataBits(9))+str2num(DataBits(11)),2);
    p2 = mod(str2num(DataBits(1))+str2num(DataBits(3))+str2num(DataBits(4))+str2num(DataBits(6))+str2num(DataBits(7))+str2num(DataBits(10))+str2num(DataBits(11)),2);
    p4 = mod(str2num(DataBits(2))+str2num(DataBits(3))+str2num(DataBits(4))+str2num(DataBits(8))+str2num(DataBits(9))+str2num(DataBits(10))+str2num(DataBits(11)),2);
    p8 = mod(str2num(DataBits(5))+str2num(DataBits(6))+str2num(DataBits(7))+str2num(DataBits(8))+str2num(DataBits(9))+str2num(DataBits(10))+str2num(DataBits(11)),2);
    p16 = mod(str2num(DataBits(1))+str2num(DataBits(2))+str2num(DataBits(3))+str2num(DataBits(4))...
        +str2num(DataBits(5))+str2num(DataBits(6))+str2num(DataBits(7))+str2num(DataBits(8))...
        +str2num(DataBits(9))+str2num(DataBits(10))+str2num(DataBits(11))+p1+p2+p4+p8, 2);
    Code = [num2str(p1) num2str(p2) DataBits(1) num2str(p4) DataBits(2:4) num2str(p8) DataBits(5:11) num2str(p16)];    
    if length(findstr(Code, '1'))==NumOnes
        ii = ii+1;
        GoodCode{ii} = Code;
    end
end

yy = GoodCode;


