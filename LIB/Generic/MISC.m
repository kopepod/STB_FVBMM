classdef MISC
 
 methods(Static = true)
  
  function dockStyle()
   set(0, 'ShowHiddenHandles', 'on');
   set(gcf,'WindowStyle','docked');
   set(0,'DefaultFigureWindowStyle','docked')
   %warning('off','all');
  end
  
  function files = filesDB(dbpath)
   ppath = pwd;
   cd(dbpath);
   filedir = struct2cell(dir);
   files = cell(1e+03,1);
   n = 0;
   for k = 1:size(filedir,2)
    if isempty(strfind(filedir{1,k},'.'))
     n = n+1;
     files{n} = filedir{1,k};
    end
   end
   files = files(1:n);
   for k = 1:n
    cd(files{k})
    x = struct2cell(dir(fullfile('*.avi')));
    x = x(1,:)';
    files{k} = x;
    cd ..
   end
   cd(ppath);
  end
    
  function bwd2uint64(dbprfx,prefix)
   for k = 1:6
    c = cell(1);
    load(strcat(dbprfx,num2str(k),'.mat'));
    for j = 1:numel(c)
     c{j}{2} = MISC.bin2uint64(c{j}{2});
    end
    save(strcat(prefix,num2str(k),'.mat'),'c','-v7.3')
   end
   
  end
  
  function y = bin2uint64(x)
   y = cell(size(x,1),1);
   for k = 1:numel(y)
    y{k} = MISC.logic2uint64(x(k,4:end));
   end
   y = cell2mat(y);
  end
  
  function x = logic2uint64(x)
   x = logical(x);
   p = uint64(2.^(63:-1:0));
   c = cell(1,size(x,2)/64);
   idx = 1;
   for j = 1:numel(c)
    s = uint64(0);
    pa = p(x(idx:idx+63));
    for k = 1:numel(pa)
     s = s + pa(k);
    end
    idx = idx+64;
    c{j} = s;
   end
   x = cell2mat(c);
  end
  
  
 end
end