classdef DTBWD
 
 properties(Constant = true)
  %{
   bn = cos(0:2*pi/6:2*pi-2*pi/6) + 1j * sin(0:2*pi/6:2*pi-2*pi/6);
   bs = {...
    logical([0,0,1]),...
    logical([0,1,0]),...
    logical([0,1,1]),...
    logical([1,0,0]),...
    logical([1,0,1]),...
    logical([1,1,0])};
  %}
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
   fields = DTBWD.getRegions([32,15]);
   o = cos([0,pi/2,pi,3*pi/2]) + 1j * sin([0,pi/2,pi,3*pi/2]);
 end
 
 methods(Static = true)
  
 
  function code = encode(v,stip)
   
   t = stip(1)-14:stip(1);
   stip = stip(2:end) - 15;% KTH stip(11:end) - 15;
   %stip = stip(11:end);
   %stip = stip(2:end);
   
   %if t(1) == 1;%WARNING REMOVED 09/FEB/2016
   % return
   %end
   
   sv = zeros(32,32,15);
   
   x = stip(1:2:end);
   y = stip(2:2:end);
   
   
   for k = 1:15
    sv(:,:,k) = v(y(k):y(k)+31,x(k):x(k)+31,t(k));
   end
   
   scan = DTBWD.fields{1,1};
   
   %sv = abs(diff(sv,[],3));
   
   sv = DTBWD.flipvideo(sv);
   
   sv = double(sv(:));
   
   tic
   %code = cellfun(@(y,x) DTBWD.subenc(sv,y,x),scan(:,1),scan(:,2));
   
   codex = cellfun(@(y,x) DTBWD.subenc(real(sv),y,x),scan(:,1),scan(:,2),'UniformOutput',false);
   codey = cellfun(@(y,x) DTBWD.subenc(imag(sv),y,x),scan(:,1),scan(:,2),'UniformOutput',false);

   codex = cell2mat(codex)';
   codey = cell2mat(codey)';   
   toc
   
   scan = scan';
   
   A = double(cell2mat(scan(1,:)));
   B = double(cell2mat(scan(2,:)));
   
   
   tic
   
   cx = mexbwd(real(sv),A,B);
   cy = mexbwd(imag(sv),A,B);
   
   toc
   
   %code = code';
   
   %x0 = x(1);
   %y0 = y(1);
   
   codexy = [false(1,22),DTBWD.BDT(x,y)];
   
   %code  = [x0,y0,t(1),codexy,code];
   
   code = [codexy,codex,codey];
   
   code = MISC.logic2uint64(code);
   
  end
  
   
  function code = BDT(x,y)
   
   x = diff(x);
   y = diff(y);
   
   %x = x/(sum(x)+eps);
   %y = y/(sum(y)+eps);

   s = sqrt(x.^2 + y.^2)+eps;
   
   u = x./s + 1j * y./s;
   
   s = s/(max(s)+eps);
   
   code = cell(1,numel(s));
   
   for k = 1:numel(s)
    
    if s(k) > 0.2 && s(k) < 0.8
     [~,idx] = min(bsxfun(@minus,u(k),DTBWD.bn));
     code{k} = DTBWD.bs{idx};
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
  
  function code = B0DT(x,y)
   %90.74 accuracy
   x = diff(x);
   y = diff(y);
   
   %x = x/(sum(x)+eps);
   %y = y/(sum(y)+eps);

   s = sqrt(x.^2 + y.^2)+eps;
   
   u = x./s + 1j * y./s;
   
   s = s/(max(s)+eps);
   
   code = cell(1,numel(s));
   
   for k = 1:numel(s)
    
    if s(k) > 0.2 && s(k) < 0.8
     [~,idx] = min(bsxfun(@minus,u(k),DTBWD.bn));
     code{k} = DTBWD.bs{idx};
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
  
  function code = encodeDTnDTBWD(v,fields,stip)
   
   t = stip(1)-14:stip(1);
   stip = round(stip(11:end) - 15);
   %stip = stip(11:end);
   
   if t(1) == 1;
    return
   end
   
   sv = zeros(32,32,15);
   
   x = stip(1:2:end);
   y = stip(2:2:end);
   
   
   
   for k = 1:15
    sv(:,:,k) = v(y(k):y(k)+31,x(k):x(k)+31,t(k));
   end
   
   scan = fields{1,1};
   
   %sv = abs(diff(sv,[],3));
   
   sv = DTBWD.flipvideo(sv);
   
   sv = double(sv(:));
   
   %tic
   code = cellfun(@(y,x) DTBWD.subenc(sv,y,x),scan(:,1),scan(:,2));

   %toc
   
   code = code';
   
   x0 = x(1);
   y0 = y(1);
   
   x(end+1) = x(end);
   y(end+1) = y(end);
   
   x = diff(x);
   y = diff(y);
   
   x = x/(sum(x)+eps);
   y = y/(sum(y)+eps);
   
   code  = [x0,y0,t(1),x,y,code];
   
   
  end
  
  
  function code = encodeDTBWD(v,fields,stip)
   
   
   s_ = floor(9*sqrt(stip(4)));
   t_ = floor(9*sqrt(stip(5)));%WARNING!!!!!
   
   scan = fields{s_,t_};
   
   y0 = stip(1);
   x0 = stip(2);
   t0 = stip(3);
   
   ry = floor((1:s_)+y0-s_/2);
   rx = floor((1:s_)+x0-s_/2);
   rt = floor((1:t_+1)+t0-t_/2);
   
   sv = v(ry,rx,rt);
   
   sv = abs(diff(sv,[],3));
   
   sv = DTBWD.flipvideo(sv);
   
   sv = double(sv(:));
   
   
   %tic
   code = cellfun(@(y,x) DTBWD.subenc(sv,y,x),scan(:,1),scan(:,2));
   
   
   %toc
   
   code = code';
   
  end
  
  function test = subenc(sv,scanx,scany)
   scanx = mexmean(sv(scanx));
   scany = mexmean(sv(scany));
   test = scanx > scany;
  end

  
  %% Patterns
  
  function regions = getRegions(st)
   
    s_ = st(1);
    t_ = st(2);
    
    v = zeros(s_,s_,t_);
    
    Y = DTBWD.splitInt(size(v,1),4);
    X = DTBWD.splitInt(size(v,2),4);
    T = DTBWD.splitInt(size(v,3),2);
    
    v = mat2cell(v,Y,X,T);
    
    v = DTBWD.fillvolumeidx(v);
    
    pttn = DTBWD.pttnfield();
    scan = cell(numel(pttn),2);
    
    n = 0;
    
    for k = 1:size(pttn,1)
     A = [pttn{k,1};pttn{k,1}+16];
     B = [pttn{k,2};pttn{k,2}+16];
     n = n+1;
     x = ismember(v,A);
     y = ismember(v,B);
     %[x,y] = BWD.eqweight(x,y);
     scan(n,:) = {x,y};
     %fprintf('%d:%f %f\n',n,sum(x)/numel(x),sum(y)/numel(y));
    end
    
    for k = 1:size(pttn,1)
     A = [pttn{k,1};pttn{k,2}+16];
     B = [pttn{k,2};pttn{k,1}+16];
     n = n+1;
     x = ismember(v,A);
     y = ismember(v,B);
     %[x,y] = BWD.eqweight(x,y);
     scan(n,:) = {x,y};
     %fprintf('%d:%f %f\n',n,sum(x)/numel(x),sum(y)/numel(y));
    end
    
    regions{1,1} = scan;
   
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
 
  
  function fields = getPatternsv2(s_,t_)
   
   pttn = DTBWD.pttnfield();
   
   %s_ = 32;
   %t_ = 15;
   
   v = zeros(s_,s_,t_);
   
   Y = DTBWD.splitInt(size(v,1),4);
   X = DTBWD.splitInt(size(v,2),4);
   T = DTBWD.splitInt(size(v,3),2);
   
   v = mat2cell(v,Y,X,T);
   
   pn = true;
   
   for t = 1:2
    c = 0;
    for x = 1:4
     for y = 1:4
      if pn
       c = c+1;
      else
       c = c-1;
      end
      sv = v{y,x,t};
      sv(:,:,:) = c;
      v{y,x,t} = sv;
     end
    end
    pn = false;
   end
   
   v = cell2mat(v);
   
   v = v(:);
   
   scan = cell(64,2);
   n = 0;
   
   for k = 1:size(pttn,1)
    
    spttn = pttn{k,1}';
    
    idx1 = false(numel(v),1);
    
    for j = spttn
     idx1 = or(idx1,abs(v) == j);
    end
    
    spttn = pttn{k,2}';
    
    idx2 = false(numel(v),1);
    
    for j = spttn
     idx2 = or(idx2,abs(v) == j);
    end
    
    n = n+1;
    
    scan(n,:) = {idx1,idx2};
    
   end
   
   
   for k = 1:size(pttn,1)
    
    idx1 = false(numel(v),1);
    
    spttn = pttn{k,1}';
    
    for j = spttn
     idx1 = or(idx1,v == j);
    end
    
    spttn = pttn{k,2}';
    
    for j = spttn
     idx1 = or(idx1,v == -j);
    end
    
    
    idx2 = false(numel(v),1);
    
    spttn = pttn{k,2}';
    
    for j = spttn
     idx2 = or(idx2,v == j);
    end
    
    spttn = pttn{k,1}';
    
    for j = spttn
     idx2 = or(idx2,v == -j);
    end
    
    
    n = n+1;
    
    scan(n,:) = {idx1,idx2};
    
   end
   
    
   fields = {scan};
   
   
   
  end
  
  
  function fields = getPatterns(stip)
   fields = cell(32);
   stip = sqrt(stip);
   
   pttn = DTBWD.pttnfield();
   
   while(~isempty(stip))
    s_ = floor(9 * stip(1,4));
    t_ = floor(9 * stip(1,5));%WARNING!!!!!
    
    v = zeros(s_,s_,t_);
    
    Y = DTBWD.splitInt(size(v,1),4);
    X = DTBWD.splitInt(size(v,2),4);
    T = DTBWD.splitInt(size(v,3),2);
    
    v = mat2cell(v,Y,X,T);
    
    pn = true;
    
    for t = 1:2
     c = 0;
     for x = 1:4
      for y = 1:4
       if pn
        c = c+1;
       else
        c = c-1;
       end
       sv = v{y,x,t};
       sv(:,:,:) = c;
       v{y,x,t} = sv;
      end
     end
     pn = false;
    end
    
    v = cell2mat(v);
    
    v = v(:);
    
    scan = cell(64,2);
    n = 0;
    
    for k = 1:size(pttn,1)
     
     spttn = pttn{k,1}';
     
     idx1 = false(numel(v),1);
     
     for j = spttn
      idx1 = or(idx1,abs(v) == j);
     end
     
     spttn = pttn{k,2}';
     
     idx2 = false(numel(v),1);
     
     for j = spttn
      idx2 = or(idx2,abs(v) == j);
     end
     
     n = n+1;
     
     scan(n,:) = {idx1,idx2};
     
    end
    
    
    for k = 1:size(pttn,1)
     
     idx1 = false(numel(v),1);
     
     spttn = pttn{k,1}';
     
     for j = spttn
      idx1 = or(idx1,v == j);
     end
     
     spttn = pttn{k,2}';
     
     for j = spttn
      idx1 = or(idx1,v == -j);
     end
     
     
     idx2 = false(numel(v),1);
     
     spttn = pttn{k,2}';
     
     for j = spttn
      idx2 = or(idx2,v == j);
     end
     
     spttn = pttn{k,1}';
     
     for j = spttn
      idx2 = or(idx2,v == -j);
     end
     
     
     n = n+1;
     
     scan(n,:) = {idx1,idx2};
     
    end    
  
    
    fields{s_,t_} = scan;
    stip(and(stip(:,4) == stip(1,4),stip(:,5) == stip(1,5)),:) = [];
   end
  end
  
  function ks = pttnscan(sv,idx)
   ks = sv(idx);
   ks = vertcat(ks{:});
  end
  
  function px = mergepx(c)
   px = c(:);
  end
  
  function vr = flipvideo(v)
   vr = v;
   v = abs(v);
   [Y,X,T] = size(v);
   hv = v(:,:,1:floor(T/2));
   [r,~] = DTBWD.xytSpV(Y,X,1);
   hv = sum(hv,3);
   hv = hv(:);
   x = hv .* r(:,1);
   y = hv .* r(:,2);
   c = sum([y,x])/(sum(hv)) - [Y,X]/2;
   c = c(1)+1j*c(2);
   c = c/abs(c);
   c = DTBWD.o-c;
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
  
  
  function sv = subVidSTIP(v,s)
   r(1) = floor(5*s(4));
   r(2) = floor(5*s(4));
   r(3) = floor(5*s(5));
   sv = v(...
    s(1)-r(1):s(1)+r(1)+mod(r(1),2)-1,...
    s(2)-r(2):s(2)+r(2)+mod(r(2),2)-1,...
    s(3)-r(3):s(3)+r(3)+mod(r(3),2)-1);
  end
  
   function field = pttnfield()
   
   field = cell(1,2);
   
   field(01,:) = {[1;2;3;4;5;6;7;8],setdiff((1:16)',[1;2;3;4;5;6;7;8])};
   field(02,:) = {[1;5;9;13;2;6;10;14],setdiff((1:16)',[1;5;9;13;2;6;10;14])};
   field(03,:) = {[1;5;2;6;11;15;12;16],setdiff((1:16)',[1;5;2;6;11;15;12;16])};
   field(04,:) = {[1;5;9;13;4;8;12;16],setdiff((1:16)',[1;5;9;13;4;8;12;16])};
   field(05,:) = {[1;2;3;4;9;10;11;12],setdiff((1:16)',[1;2;3;4;9;10;11;12])};
   
   field(06,:) = {[1;2;3;4;13;14;15;16],setdiff((1:16)',[1;2;3;4;13;14;15;16])};
   field(07,:) = {[1;5;9;13;3;7;11;15],setdiff((1:16)',[1;5;9;13;3;7;11;15])};
   field(08,:) = {[2;6;3;7;9;13;12;16],setdiff((1:16)',[2;6;3;7;9;13;12;16])};
   field(09,:) = {[1;2;7;8;11;12;13;14],setdiff((1:16)',[1;2;7;8;11;12;13;14])};
   field(10,:) = {[1;4;6;7;10;11;13;16],setdiff((1:16)',[1;4;6;7;10;11;13;16])};
   
   
   field(11,:) = {[1;5;3;7;10;14;12;16],setdiff((1:16)',[1;5;3;7;10;14;12;16])};%corrected
   field(12,:) = {[1;2;7;8;9;10;15;16],setdiff((1:16)',[1;2;7;8;9;10;15;16])};
   field(13,:) = {[1;3;6;10;8;12;13;15],setdiff((1:16)',[1;3;6;10;8;12;13;15])};%corrected
   field(14,:) = {[1;4;6;7;9;12;14;15],setdiff((1:16)',[1;4;6;7;9;12;14;15])};%corrected
   field(15,:) = {[1;3;6;8;9;11;14;16],setdiff((1:16)',[1;3;6;8;9;11;14;16])};

%    field(16,:) = {[1;2;3;4;5;8;10;11],setdiff((1:16)',[1;2;3;4;5;8;10;11])};%WARNING!!!!!
%    field(17,:) = {[1;2;5;7;9;11;13;14],setdiff((1:16)',[1;2;5;7;9;11;13;14])};
% 
%    field(30,:) = {[1;2;3;4;6;7],[10;11;13;14;15;16]};
%    field(31,:) = {[1;5;6;9;10;13],[4;7;8;11;12;16]};
%    field(32,:) = {[1;2;3;4;5;8],[9;12;13;14;15;16]};
%    field(33,:) = {[2;3;5;6;7;8],[9;10;11;12;14;15]};
%    field(34,:) = {[2;5;6;9;10;14],[3;7;8;11;12;15]};
%    
%    field(35,:) = {[1;2;3;5;6;9],[8;11;12;14;15;16]};
%    field(36,:) = {[2;3;4;7;8;12],[5;9;10;13;14;15]};
%    field(37,:) = {[1;2;3;4;10;11],[6;7;13;14;15;16]};
%    field(38,:) = {[1;4;5;6;7;8],[9;10;11;12;13;16]};
%    field(39,:) = {[1;2;5;9;13;14],[3;4;8;12;15;16]};
%    
%    field(40,:) = {[5;6;7;8;14;15],[2;3;9;10;11;12]};
%    field(41,:) = {[1;5;7;9;11;13],[4;6;8;10;12;16]};
%    field(42,:) = {[1;2;6;10;13;14],[3;4;7;11;15;16]};
%    field(43,:) = {[2;6;8;10;12;14],[3;5;7;9;11;15]};   

   field(60,:) = {[1;2;5;6],[3;7;4;8]};
   field(61,:) = {[9;13;10;14],[11;15;12;16]};
   field(62,:) = {[1;2;3;4],[5;6;7;8]};
   field(63,:) = {[9;10;11;12],[13;14;15;16]};
   field(64,:) = {[1;5;3;7],[2;6;4;8]};

   field(65,:) = {[10;14;12;16],[9;13;11;15]};
   field(66,:) = {[1;5;2;6],[9;13;10;14]};
   field(67,:) = {[3;7;4;8],[11;15;12;16]};
   field(68,:) = {[1;5;9;13],[2;6;10;14]};
   field(69,:) = {[3;7;11;15],[4;8;12;16]};
 
   field(70,:) = {[1;2;9;10],[5;6;13;14]};
   field(71,:) = {[3;4;11;12],[7;8;15;16]};
   field(72,:) = {[1;5;2;6],[11;15;12;16]};
   field(73,:) = {[3;7;4;8],[9;10;13;14]};
   field(74,:) = {[1;5;4;8],[10;11;14;15]};

   field(75,:) = {[1;2;13;14],[7;11;8;12]};
   field(76,:) = {[2;3;14;15],[5;9;8;12]};
   
   
%    field(77,:) = {[5;6;9;10],[7;8;10;12]};%-added
%    field(78,:) = {[5;6;7;8],[9;10;11;12]};
%    field(79,:) = {[2;3;6;7],[10;11;14;15]};
%    field(80,:) = {[2;6;10;14],[3;7;11;15]};


   
   for k = size(field,1):-1:1
    if isempty(field{k,1})
     field(k,:) = [];
    end
   end
   
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
  
 end
end