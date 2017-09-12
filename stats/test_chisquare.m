% Test for conditional indepndence via chi-square test (from Kevin Murphy's BNT)
%
% test CHI2
%

% Génération de données jouet
N=20000;
xx=5*rand(N,1); 
yy=5*rand(N,1); 
zz = 5*rand(N,1);
z2= (xx+yy) ; 
Data=1+[round(xx) round(yy) round(z2) round(zz)];

texte{1}='dependency';
texte{2}='independency';

fprintf('\n\n\n\t=================================================\n');
fprintf('\t\tChi2 Test for (conditional) independency\n');
fprintf('\t\t(X, Y and Z are theoritically independant)\n');

fprintf('\t=================================================\n');
fprintf('Variables \t\t: Theory :   Results\n');
fprintf('------------------------------------------------------------------------\n');

S=[];
x=1; y=2;
[IND CHI2] = cond_indep_chisquare(x,y,S,Data(1:100,:));
fprintf('X, Y \t\t\t: DEPDCY : \tCHI2=%8.2f \t%s (not enough data)\n',CHI2,texte{IND+1});

S=[];
x=1; y=2;
[IND CHI2] = cond_indep_chisquare(x,y,S,Data);
fprintf('X, Y \t\t\t: INDPCY : \tCHI2=%8.2f \t%s\n',CHI2,texte{IND+1});

S=[];
x=1; y=3;
[IND CHI2] = cond_indep_chisquare(x,y,S,Data);
fprintf('X, X+Y \t\t\t: DEPDCY :\tCHI2=%8.2f \t%s\n',CHI2,texte{IND+1});

S=[];
x=1; y=4;
[IND CHI2] = cond_indep_chisquare(x,y,S,Data);
fprintf('X, Z \t\t\t: INDPCY :\tCHI2=%8.2f \t%s\n',CHI2,texte{IND+1});


x=1; y=2; S=4;
[IND CHI2] = cond_indep_chisquare(x,y,S,Data);
fprintf('X, Y | Z \t\t: INDPCY : \tCHI2=%8.2f \tconditional %s\n',CHI2,texte{IND+1});

x=1; y=2; S=3;
[IND CHI2] = cond_indep_chisquare(x,y,S,Data);
fprintf('X, Y | X+Y \t\t: DEPDCY : \tCHI2=%8.2f \tconditional %s\n',CHI2,texte{IND+1});

x=1; y=2; S=[3 4];
[IND CHI2] = cond_indep_chisquare(x,y,S,Data);
fprintf('X, Y | {X+Y,Z} \t: DEPDCY : \tCHI2=%8.2f \tconditional %s\n',CHI2,texte{IND+1});

fprintf('------------------------------------------------------------------------\n');