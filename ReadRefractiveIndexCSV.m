function [n,k,wls]=ReadRefractiveIndexCSV(name,R1,C1,R2,C2)
%read n,k data from refractiveindex.info
%format of csv must be like Si.csv
%[wls,n] delimiter: row(R1:R2),colum(C1:C2)
A=csvread(name,R1-1,C1-1,[R1-1,C1-1,R2-1,C2-1]);
n=A(:,2);
wls=A(:,1);
B=csvread(name,R1+R2,C1-1,[R1+R2,C1-1,R2+R2,C2-1]);
k=B(:,2);
end