function avgHrfGamma(loadfile)
    load(loadfile,'T','W','A','r','r2','hemoPred')
    
    
    y = hrfGamma(t,T,W,A)

end