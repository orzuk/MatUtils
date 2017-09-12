function x=mgd(N,d,rmean,covariance)
%Generates a Multivariate Gaussian Distribution
%Usage x=mgd(N,d,mu,sigmax)
%
%This function generates N samples from a d-dimension
%Gaussian distribution. It differs from randn by the fact that
%you can specify the mean and the covariance you want when the samples
%are generated. 
%
%N - Number of samples to generate
%d - dimension of the Gaussian distribution
%mu - mean around which to center the samples;
%sigmax - covariance matrix for the samples, should be positive definite
%
%Ex. Generate 50 samples from a 2 dimensional multivariate Gaussian 
%distribution, with a mean of [4,5] and a covariance matrix of [9,0;0,9].
%x=mgd(50,2,[4,5]',[9,0;0,9]) or x=mgd(50,2,[4 5],[9 0;0 9])
%new version of mgd doesn't use whitening transform
%Updated May 17 2007, Updated Feb 5th 2006, original version 2004.
%By Timothy Felty - Bug Fix By Or Zuk

%Parse input parameters
[rowsm,colsm]=size(rmean);
lengthmean=length(rmean);
[rowsv,colsv]=size(covariance);
if rowsv ~= colsv 
    error('Covariance matrix should be square')
end
if lengthmean ~= rowsv
    error('The dimension of the requested mean is not equal to the dimensions of the covariance')
end
if d ~= lengthmean
    error('The mean should have the same dimension as the requested samples')
end
if N < 1
    N=1;
end
N=fix(N); %Make sure N is a whole number

%this allows for generality when passing the mean into the function
if (colsm==1)  %if its a column vector turn it into a row vector
    rmean=rmean';
end

%Start of the program
x=randn(N,d);  %generate the samples using built in Matlab function

%This makes the Matlab distribution have the covariance wanted by the user
%Computes the Cholesky decomposition of the given variance
%It uses a different method depending on wheter or not variance is positive
%definite. This is used because the the variance is needed which is the
%square root of the covariance
[R,p]=chol(covariance);
if p>0
    x=x*sqrtm(covariance);
else
    x=x*R;
end

%Add back on the desired mean
x=x+repmat(rmean,[N,1]);

    
