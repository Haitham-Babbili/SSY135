function mod= get_modulation(M)
%%

switch M
    case 4
      mod =[1+1i 1-1i -1-1i -1+1i]/sqrt(2); 
        
    otherwise
        error('ERROR in get_modulation, there is no implemented modulation for this M')
 

end