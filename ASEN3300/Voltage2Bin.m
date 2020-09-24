function [binnumout,binaryout]=Voltage2Bin(min_voltage, max_voltage, bits,voltage)

range = max_voltage - min_voltage;
N = 2^(bits);

LSB = range/N;

binnumout = floor((voltage-min_voltage)./LSB);
binaryout = de2bi(abs(binnumout));
end