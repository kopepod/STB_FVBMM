
for n = [100,200,400,600,800,1200,1500,2000,4000]

 for C = [1e-3,1e-2,1e-1,1,1e1,1e2,1e3]
  
  [labels, A] = nAkELM(TX', TY', LX, C, n);
  
  disp('__________')
  disp([n,C])
  disp(sum(labels == LY) / numel(LY))
  
  
 end
 
end


