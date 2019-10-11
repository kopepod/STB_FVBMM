classdef BDTMBH
 
 properties(Constant = true)
  bn = cos(0:2*pi/12:2*pi-2*pi/12) + 1j * sin(0:2*pi/12:2*pi-2*pi/12);
   bs = {...
    logical([0,0,1]),...
    logical([0,1,0]),...
    logical([0,1,1]),...
    logical([1,0,0]),...
    logical([1,0,1]),...
    logical([1,1,0]),...
    logical([0,0,1]),...
    logical([0,1,0]),...
    logical([0,1,1]),...
    logical([1,0,0]),...
    logical([1,0,1]),...
    logical([1,1,0])};
   o = cos([0,pi/2,pi,3*pi/2]) + 1j * sin([0,pi/2,pi,3*pi/2]);
   fields = BDTMBH.getRegions([32,15]);
   yx = BDTMBH.xytSpV(32,32,1);
 end
 
 methods(Static = true)
  
 
  function code = encode(v,stip)
   
   t = stip(1)-14:stip(1);
   stip = stip(2:end) - 15;
   
   [~,~,T] = size(v);
   
   sv = zeros(32,32,15);
   
   x = stip(1:2:end);
   y = stip(2:2:end);
   
   t(t < 1) = 1;
   t(t > T) = T;
   
   for k = 1:15
    ry = y(k):y(k)+31;
    rx = x(k):x(k)+31;
    sv(:,:,k) = v(ry,rx,t(k));
   end
   
   scan = BDTMBH.fields;
   
   svx = BDTMBH.flipvideo(real(sv));
   svy = BDTMBH.flipvideo(imag(sv));
   
   svx = svx(:);
   svy = svy(:);
   
   codex = mexbwd(svx,scan{1},scan{2});
   codey = mexbwd(svy,scan{1},scan{2});
   
   T = BDTMBH.binTraj(x,y);
   T = MISC.logic2uint64([false(1,22),T]);
   
   code  = [T,codex,codey];
   
  end
   
  function code = binTraj(x,y)
   
   x = diff(double(x));
   y = diff(double(y));
   
   s = sqrt(x.^2 + y.^2)+eps;
   
   u = x./s + 1j * y./s;
   
   s = s/(max(s)+eps);
   
   code = cell(1,numel(s));
   
   for k = 1:numel(s)
    
    if s(k) > 0.2 && s(k) < 0.8
     [~,idx] = min(bsxfun(@minus,u(k),BDTMBH.bn));
     code{k} = BDTMBH.bs{idx};
     continue
    end
    
    if s(k) > 0.2
     code{k} = true(1,3);
    else
     code{k} = false(1,3);
    end
    
    
   end
   
   code = cell2mat(code);
   
  end

  function vr = flipvideo(v)
   vr = v;
   [Y,X,T] = size(v);
   hv = v(:,:,1:floor(T/2));
   hv = sum(hv,3);
   hv = hv(:);
   x = hv .* BDTMBH.yx(:,1);
   y = hv .* BDTMBH.yx(:,2);
   c = sum([y,x])/(sum(hv)) - [Y,X]/2;
   c = c(1)+1j*c(2);
   c = c/abs(c);
   c = BDTMBH.o-c;
   c = c .* conj(c);
   [~,c] = min(c);
   switch c
    case 1
     vr =  fliplr(vr);
    case 2
     vr =  flipud(vr);
    otherwise
    return
   end
  end
 
  function [r,p] = xytSpV(nx,ny,nt)
   p = struct;
   p.x = zeros(nx,ny,nt);
   p.y = zeros(nx,ny,nt);
   p.t = zeros(nx,ny,nt);
   for t = 1:nt
    for x = 1:nx
     for y = 1:ny
      p.y(y,x,t) = y;
      p.x(y,x,t) = x;
      p.t(y,x,t) = t;
     end
    end
   end
   a = p.y;
   b = p.x;
   c = p.t;
   r = [a(:),b(:),c(:)];
  end
 
  %% Scanning Regions
  
  function regions = getRegions(st)
   
   s_ = st(1);
   t_ = st(2);
   
   v = zeros(s_,s_,t_);
   
   Y = BDTMBH.splitInt(size(v,1),4);
   X = BDTMBH.splitInt(size(v,2),4);
   T = BDTMBH.splitInt(size(v,3),2);
   
   v = mat2cell(v,Y,X,T);
   
   v = BDTMBH.fillvolumeidx(v);
   
   pttn = BDTMBH.pttnfield();
   scan = cell(1e2,2);
   
   n = 0;
   
   %space symmetric
   for k = 1:32%size(pttn,1)
    A = [pttn{k,1};pttn{k,1}+16];
    B = [pttn{k,2};pttn{k,2}+16];
    n = n+1;
    x = ismember(v,A);
    y = ismember(v,B);
    scan(n,:) = {x,y};
   end
   
   % time symmetric
   %%{
   for k = 1:32%size(pttn,1)
    A = [pttn{k,1};pttn{k,2}+16];
    B = [pttn{k,2};pttn{k,1}+16];
    n = n+1;
    x = ismember(v,A);
    y = ismember(v,B);
    scan(n,:) = {x,y};
   end
   
   scan = scan(1:n,:)';
   
   scan = {cell2mat(scan(1,:)),cell2mat(scan(2,:))};
   
   regions = scan;
  
  end
  
  function field = pttnfield()
   
   field = cell(32,2);
   
   field(01,:) = {[1;2;3;4;5;6;7;8],[9;10;11;12;13;14;15;16]};
   field(02,:) = {[1;5;9;13;2;6;10;14],[3;7;11;15;4;8;12;16]};
   field(03,:) = {[1;5;2;6;11;15;12;16],[9;13;10;14;3;7;4;8]};
   field(04,:) = {[1;5;9;13;4;8;12;16],[2;6;10;14;3;7;11;15]};
   field(05,:) = {[1;2;3;4;9;10;11;12],[5;6;7;8;13;14;15;16]};
   
   field(06,:) = {[1;2;3;4;13;14;15;16],[5;6;7;8;9;10;11;12]};
   field(07,:) = {[1;5;9;13;3;7;11;15],[2;6;10;14;4;8;12;16]};
   field(08,:) = {[2;6;3;7;9;13;12;16],[1;5;10;11;14;15;4;8]};
   field(09,:) = {[1;2;7;8;11;12;13;14],[3;4;5;6;9;10;15;16]};
   field(10,:) = {[1;4;6;7;10;11;13;16],[2;3;5;9;8;12;14;15]};
   
   
   field(11,:) = {[1;5;3;7;13;14;12;16],[2;6;4;8;9;13;11;15]};
   field(12,:) = {[1;2;7;8;9;10;15;16],[3;4;5;6;11;12;13;14]};
   field(13,:) = {[1;3;6;10;8;12;13;15],[2;4;5;6;7;11;14;16]};
   field(14,:) = {[1;4;6;7;9;11;14;16],[2;4;5;7;10;12;13;15]};
   field(15,:) = {[1;3;6;8;9;11;14;16],[2;4;5;7;10;12;13;15]};
   
   field(16,:) = {[1;2;5;6],[3;7;4;8]};
   field(17,:) = {[9;13;10;14],[11;15;12;16]};
   field(18,:) = {[1;2;3;4],[5;6;7;8]};
   field(19,:) = {[9;10;11;12],[13;14;15;16]};
   field(20,:) = {[1;5;3;7],[2;6;4;8]};
   
   field(21,:) = {[10;14;12;16],[9;13;11;15]};
   field(22,:) = {[1;5;2;6],[9;13;10;14]};
   field(23,:) = {[3;7;4;8],[11;15;12;16]};
   field(24,:) = {[1;5;9;13],[2;6;10;14]};
   field(25,:) = {[3;7;11;15],[4;8;12;16]};
   
   field(26,:) = {[1;2;9;10],[5;6;13;14]};
   field(27,:) = {[3;4;11;12],[7;8;15;16]};
   field(28,:) = {[1;5;2;6],[11;15;12;16]};
   field(29,:) = {[3;7;4;8],[9;10;13;14]};
   field(30,:) = {[1;5;4;8],[10;11;14;15]};
   
   field(31,:) = {[1;2;13;14],[7;11;8;12]};
   field(32,:) = {[2;3;14;15],[5;9;8;12]};
   
   
  end
  
  function v = fillvolumeidx(v)
   idx = 0;
   for t = 1:size(v,3)
    for x = 1:size(v,2)
     for y = 1:size(v,1)
      idx = idx+1;
      sv = v{y,x,t};
      sv(:,:,:) = idx;
      v{y,x,t} = sv;
     end
    end
   end
   v = cell2mat(v);
   v = v(:);
  end

  function s = splitInt(a,n)
   x = floor(a/n);
   s = x * ones(n,1);
   m = 0;
   while sum(s) ~= a
    m = m+1;
    s(m) = s(m)+1;
   end
  end
 
  %% video gradient source
  
  function v = videoGsrc(v)
   v = abs(diff(v,[],3));
   v(:,:,end+1) = v(:,:,end);
   v = double(v);
   vx = convn(v,[-1,0,1],'same');
   vy = convn(v,[-1;0;1],'same');
   v = abs(vx)+1j*abs(vy);
  end
  
 end
end