function x = UDPdecon(UDP, UDPgluc)

x0 = [0.5 0 0 0 0 0 0.5]; %starting guess

[x, resnorm] = lsqcurvefit(@MIDfit,x0,UDP,UDPgluc, [0 0 0 0 0 0 0],[1 1 1 1 1 1 1])%least squares fit to find glucose labelling in UDP glucose

    function cc = MIDfit(b,a)
        
        
        %correction matricies which map combination of two compounds into
        %one labellng pattern 
        cora = [2,1,1,1,1,1,1;2,3,1,1,1,1,1;2,3,4,1,1,1,1;2,3,4,5,1,1,1;2,3,4,5,6,1,1;2,3,4,5,6,7,1;2,3,4,5,6,7,8;3,4,5,6,7,8,9;4,5,2,7,8,9,10;5,6,7,8,9,10,11;6,7,8,9,10,11,1;7,8,9,10,11,1,1;8,9,10,11,1,1,1;9,10,11,1,1,1,1;10,11,1,1,1,1,1;11,1,1,1,1,1,1];
        corb = [2,1,1,1,1,1,1;3,2,1,1,1,1,1;4,3,2,1,1,1,1;5,4,3,2,1,1,1;6,5,4,3,2,1,1;7,6,5,4,3,2,1;8,7,6,5,4,3,2;8,7,6,5,4,3,2;8,7,6,5,4,3,2;8,7,6,5,4,3,2;8,7,6,5,4,3,1;8,7,6,5,4,1,1;8,7,6,5,1,1,1;8,7,6,1,1,1,1;8,7,1,1,1,1,1;8,1,1,1,1,1,1];
        
        %add extra zero to account for some empy spaces in correction
        %matrix
        a = [0 a]';
        b = [0 b]';
        
        %populate correction matricies with MIDs
        aa = a(cora);
        bb = b(corb);
        
        %multiply correction matricies to combine labelling patterns 
        c = aa.*bb;
        
        %sum matrix to get MID of combined molecule 
        cc = sum(c');
    end
end