function [CI, Chi2] = cond_indep_chisquare(X, Y, S, Data, alpha)
% COND_INDEP_CHISQUARE Test if X indep Y given Z using ChiSquare test
% [CI Chi2] = cond_indep_chisquare(X, Y, S, Data, alpha)
%
% Input : 
%  	Data is the data matrix, N rows * NbVar columns
%	X is the index of variable X in Data matrix
%	Y is the index of variable Y in Data matrix
%	S are the indexes of variables in set S
% 	alpha is the significance level (default: 0.01)
%
% Output :
%	CI = test result (1=conditional independency, 0=no)
%	Chi2 = chi2 value (-1 if not enough data to perform the test --> CI=0)
%
%
% V1.3 : 14 feb 2003 (Ph. Leray)
%
%
% Things to do :
% - do not use 'find' in nij computation (when S=empty set)
% - find a better way than 'warning off/on' in tmpij, tmpijk computation
%

if nargin < 5, alpha = 0.05; end

ns = max(Data);
N = size(Data,1);
qi=ns(S);
tmp=[1 cumprod(qi(1:end-1))];
qs=1+(qi-1)*tmp';
if isempty(qs),
    nij=zeros(ns(X),ns(Y));
    df=prod(ns([X Y])-1)*prod(ns(S));
else
    nijk=zeros(ns(X),ns(Y),qs);
    tijk=zeros(ns(X),ns(Y),qs);
    df=prod(ns([X Y])-1)*qs;
end


if (N<10*df)
    % Not enough data to perform the test
    Chi2=-1; 
    CI=0;
    
elseif isempty(S)
    for i=1:ns(X),
        for j=1:ns(Y),
            nij(i,j)=length(find((Data(:,X)==i)&(Data(:,Y)==j))) ;
        end
    end
    restr=find(sum(nij,1)==0);
    if ~isempty(restr) 
        nij=nij(:,find(sum(nij,1)));
    end
    
    tij=sum(nij,2)*sum(nij,1)/N ;
    tmpij=nij-tij;
    
    % Correction de Yates
    [xi yj]=find(tij<10) ;
    for i=1:length(xi),
        tmpij(xi(i),yj(i))=abs(tmpij(xi(i),yj(i)))-0.5;
    end

    warning off;
    tmp=(tmpij.^2)./tij;
    warning on;
    tmp(find(tmp==Inf))=0;

    Chi2=sum(sum(tmp));

    alpha2=1-chisquared_prob(Chi2,df);
    CI=(alpha2>=alpha) ;
    
else
    for exemple=1:N,
        i=Data(exemple,X);
        j=Data(exemple,Y);
        Si=Data(exemple,S)-1; 
        k=1+Si*tmp';
        nijk(i,j,k)=nijk(i,j,k)+1;
    end
    
    nik=sum(nijk,2);
    njk=sum(nijk,1);
    N2=sum(njk);
    for k=1:qs, 
        if N2(:,:,k)==0
                tijk(:,:,k)=0;
            else
                tijk(:,:,k)=nik(:,:,k)*njk(:,:,k)/N2(:,:,k);
            end
    end
    
    tmpijk=nijk-tijk;
    
    % Correction de Yates
    [xi yj]=find(tijk<10) ;
    for i=1:length(xi),
        tmpijk(xi(i),yj(i))=abs(tmpijk(xi(i),yj(i)))-0.5;
    end

    warning off;
    tmp=(tmpijk.^2)./tijk;
    warning on;
    tmp(find(tmp==Inf))=0;
    
    Chi2=sum(sum(sum(tmp)));
    alpha2=1-chisquared_prob(Chi2,df);
    CI=(alpha2>=alpha) ;
        
end