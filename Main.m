addpath(genpath(strcat(pwd,'/LIB')))
addpath(genpath(strcat(pwd,'/DATA')))

filedb = 'UCF101_BIDT_MBH';
splitmask = [];load('UCF101_splitmask');


disp(datetime)

nsplit = 2;

fprintf('split: %d\n',nsplit);

idx = splitmask{nsplit};

[lambda,F] = BMMu64.FisherBMM(filedb,idx,4950,256);

X = BMMu64.encodeBMMfisher(filedb,lambda,F);

Y = X(~idx,:);
X = X(idx,:);

LX = cell2mat(X(:,3));
X  = cell2mat(X(:,2));

LY = cell2mat(Y(:,3));
Y  = cell2mat(Y(:,2));

for M = 512
  T = pca(X,'NumComponents',M);
  
  TX = X * T;
  TY = Y * T;
  for n = 0.6
    for C = 50
      [labels, A] = nAkELM(shiftdim(TX,1), shiftdim(TY,1), shiftdim(LX,1), C, n);
      acc = sum(labels == LY) / numel(LY);
      fprintf('{%d,%0.2f,%d,%0.4f} \n',M,n,C,acc);
    end
  end
end

disp(datetime)

