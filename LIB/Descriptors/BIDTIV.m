classdef BIDTIV
 
 properties(Constant = true)
  bn = cos(0:2*pi/12:2*pi-2*pi/12) + 1j * sin(0:2*pi/12:2*pi-2*pi/12);
  bs = {...
   logical([0,0,0,1]),...
   logical([0,0,1,0]),...
   logical([0,0,1,1]),...
   logical([0,1,0,0]),...
   logical([0,1,0,1]),...
   logical([0,1,1,0]),...
   logical([0,0,0,1]),...
   logical([0,0,1,0]),...
   logical([0,0,1,1]),...
   logical([0,1,0,0]),...
   logical([0,1,0,1]),...
   logical([0,1,1,0])};
  log2u64 = uint64(2.^(63:-1:0));
  
  o = cos([0,pi/2,pi,3*pi/2]) + 1j * sin([0,pi/2,pi,3*pi/2]);
  fields = BIDTIV.getRegions([32,14]);
  yx = BIDTIV.xytSpV(32,32,1);
  visualizepairs = false;
 end
 
 methods(Static = true)
  
 
  function code = encode(v,stip)
   
   t = stip(1)-14:stip(1);
   stip = stip(2:end) - 15;
   
   [~,~,T] = size(v);
   
   sv = zeros(32,32,14);
   
   x = stip(1:2:end);
   y = stip(2:2:end);
   
   t(t < 1) = 1;
   t(t > T) = T;
   
   for k = 1:14
    ry = y(k):y(k)+31;
    rx = x(k):x(k)+31;
    sv(:,:,k) = v(ry,rx,t(k));
   end
   
   %scan = BIDTIV.fields;
   
   nregions = BIDTIV.fields;
   
   svx = BIDTIV.flipvideo(real(sv));
   svy = BIDTIV.flipvideo(imag(sv));
   
   svx = svx(:);
   svy = svy(:);
   
   code = mexBWDiv(svx,svy,nregions);
   
   
   %codex = mexbwd(svx,scan{1},scan{2});
   %codey = mexbwd(svy,scan{1},scan{2});
   
   T = [false(1,8),BIDTIV.binTraj(x',y')];
   %T = MISC.logic2uint64(T);
   
   T = sum(BIDTIV.log2u64(T),'native');
   
   
   %code  = [T,codex,codey];
   
   code  = [T,code];
   
  end
  
  function code = binTraj(x,y)
   x = diff(double(x));
   y = diff(double(y));
   s = sqrt(x.^2 + y.^2)+eps;
   u = x./s + 1j * y./s;
   s = s/(max(s));
   A = bsxfun(@minus,u,BIDTIV.bn);
   [~,idx] = min(abs(A),[],2);
   code = BIDTIV.bs(idx);
   
   idx = find(s > 0.8)';
   
   for k = idx
    code{k} = true(1,4);
   end

   idx = find(s < 0.2)';
   
   for k = idx
    code{k} = false(1,4);
   end

   code = cell2mat(code);
  end
  
  function vr = flipvideo(v)
   vr = v;
   [Y,X,T] = size(v);
   hv = v(:,:,1:floor(T/2));
   hv = sum(hv,3);
   hv = hv(:);
   x = hv .* BIDTIV.yx(:,1);
   y = hv .* BIDTIV.yx(:,2);
   c = sum([y,x])/(sum(hv)) - [Y,X]/2;
   c = c(1)+1j*c(2);
   c = c/abs(c);
   c = BIDTIV.o-c;
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
   
   Y = BIDTIV.splitInt(size(v,1),4);
   X = BIDTIV.splitInt(size(v,2),4);
   T = BIDTIV.splitInt(size(v,3),2);
   
   v = mat2cell(v,Y,X,T);
   
   v = BIDTIV.fillvolumeidx(v);
   
   pttn = BIDTIV.pttnfield();
   scan = cell(1e2,1);
   
   n = 0;
   
   c = cell(64,1);
   
   w = 0.5*ones(s_+ 2,s_+ 2,t_);
   
   %space symmetric
   for k = 1:32%size(pttn,1)
    A = [pttn{k,1};pttn{k,1}+16];
    B = [pttn{k,2};pttn{k,2}+16];
    n = n+1;
    x = ismember(v(:),A);
    y = ismember(v(:),B);
    u = zeros(size(v));
    u(:,:,:) = -1;
    u(x) = 0;
    u(y) = 1;
    w(2:s_+1,2:s_+1,:) = u;
    w(w == -1) = 0.25;
    c{n} = w;
    %{
    implay(u(1:1/8:end,1:1/8:end,:))
    %}
    scan{n} = BIDTIV.getVI(u);
   end
   
   % time symmetric
   %%{
   for k = 1:32%size(pttn,1)
    A = [pttn{k,1};pttn{k,2}+16];
    B = [pttn{k,2};pttn{k,1}+16];
    n = n+1;
    x = ismember(v(:),A);
    y = ismember(v(:),B);
    u = zeros(size(v));
    u(:,:,:) = -1;
    u(x) = 0;
    u(y) = 1;
    w(2:s_+1,2:s_+1,:) = u;
    w(w == -1) = 0.25;
    c{n} = w;
    scan{n} = BIDTIV.getVI(u);
   end
   
   
   c = reshape(c,[8,8]);
   c = cell2mat(c);
   
   if BIDTIV.visualizepairs
    implay(c)
   end
   
   scan = scan(1:n);
   
   regions = cell2mat(scan);
   
   regions = BIDTIV.mapregions4mex(regions);
  
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
   
   
   field(11,:) = {[1;5;3;7;10;14;12;16],[2;6;4;8;9;13;11;15]};
   field(12,:) = {[1;2;7;8;9;10;15;16],[3;4;5;6;11;12;13;14]};
   field(13,:) = {[1;3;6;10;8;12;13;15],[2;4;5;7;9;11;14;16]};
   field(14,:) = {[1;4;6;7;9;12;14;15],[2;3;5;8;10;11;13;16]};
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
  
  %% Explore volume

  
  function bb = getVI(u)
   bbx = cell(32,1);
   bby = cell(32,1);
   n = 0;
   m = 0;
   while numel(unique(u)) ~= 1
    [y0,x0,t0,val] = BIDTIV.findstart(u);
    [y1,x1,t1,u] = BIDTIV.boundbox(y0,x0,t0,u,val);
    if val
     n = n +1;
     bbx{n} = BIDTIV.yxt2bbox([y0,x0,t0;y1,x1,t1]);
    else
     m = m +1;
     bby{m} = BIDTIV.yxt2bbox([y0,x0,t0;y1,x1,t1]);
    end
   end
   bbx = cell2mat(bbx(1:n));
   [~,idx] = sort(bbx(:,1),'descend');
   bbx = bbx(idx,:);       
   bby = cell2mat(bby(1:m));
   [~,idx] = sort(bby(:,1),'descend');
   bby = bby(idx,:);
   bb = [bby;[0,0,0];bbx;[0,0,0]];
  end
  
  function bbox = yxt2bbox(yxt)
   bbox = zeros(8,3);
   n = 0;
   for k = 1:2
    t = yxt(k,3);
    for j = 1:2
     x = yxt(j,2);
     for  i = 1:2
      y = yxt(i,1);   
      n = n+1;
      if any([1,4,6,7] == n)
       bbox(n,:) = -[y,x,t];
      else
       bbox(n,:) = [y,x,t];
      end
     end
    end
   end
   for k = n:-1:1
    if any(abs(bbox(k,:)) == 0)
     bbox(k,:) = [];
    end
   end
  end
  
  function [y,x,t,val] = findstart(u)
   [Y,X,T] = size(u);
   for t = 0:T-1
    for x = 0:X-1
     for y = 0:Y-1
      val = u(y+1,x+1,t+1);
      if val ~= -1       
       return
      end
     end
    end
   end
   
  end
 
  function [y1,x1,t1,u] = boundbox(y0,x0,t0,u,val)
   
   [Y,X,T] = size(u);
   
   for k = 1:3    
    y = y0+1;
    x = x0+1;
    t = t0+1;
    switch k
     case 1
      while y ~= Y
       y = y+1;
       if val ~= u(y,x,t)
        y = y-1;
        break
       end
      end
      y1 = y;
     case 2
      while x ~= X
       x = x+1;
       if val ~= u(y,x,t)
        x = x-1;
        break
       end
      end
      x1 = x;
     case 3
      while t ~= T
       t = t+1;
       if val ~= u(y,x,t)
        t = t-1;
        break
       end
      end
      t1 = t;
    end
   end
   
   u(y0+1:y1,x0+1:x1,t0+1:t1) = -1;
   
  end
  
  function nregions = mapregions4mex(regions)
   a = regions(:,1);
   for k = 1:numel(a)
    if a(k) > 0
     a(k) = 1;
     regions(k,:) = abs(regions(k,:))-1;
     continue
    end
    if a(k) < 0
     a(k) = 2;
     regions(k,:) = abs(regions(k,:))-1;
     continue
    end  
   end
   nregions = uint64([a,regions]);
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