function [Du] = Function_Dose_Calculation(u)
    En = size(u,5);
    [Energy_vec] = Function_Energy_vec(En);
    Du = sum(trapz(Energy_vec,u,5),4);
end