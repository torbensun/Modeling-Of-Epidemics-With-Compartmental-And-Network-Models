 def isbetter(new_eval,current_eval):
 	if (new_eval>current_eval):
 		return 1:
 	else:
 		return 0;
 
 #must substitute f with the name of the function later 
	
 def optimize(paramenters,d_paramenters, current_itter, current_eval, new_eval, sign, isbetter,current_paramenter,last_success,O_err, og_paramenters):
 	#basecase
	if (current_itter-last_success>(O_err+6)): # failed to improve for a lot of tries
	
		if (current_paramenter<len(paramenters): #and componets are not set yet
			paramenters[current_paramenter]=paramenters[current_paramenter]+d_paramenters[current_paramenter]*((current_itter-last_success)*10**(-O_err)*sign) #undo changes to component since last success
			O_err=O_err+1 # increase precesion
			current_itter=current_itter+1
			
			if (O_err>6): 		# if your adjustment are really small and stil do nothing
				sign=-1 
				current_itter=0
				last_success=0
				O_err=1
				return f; 				# try the other direction
			 
			if (sign==-1 and O_err>5)#if you tried both directions and it doesn't improve anymore		
			d_paramenters=paramenters[current_paramenter]-og_paramenters[current_paramenter] #remeber what to change
			paramenters[current_paramenter]=og_paramenters[current_paramenter] #reset parameters
			#move to the next component
				sign=1 
				current_itter=0
				last_success=0
				O_err=0
				current_paramenter=current_paramenter+1 	
				return f;
			return f; # repeat after increasing precesion
				
		if current_paramenter>=len(paramenters): #if all components are set
		
			if (O_err>6): # and it's really minute
				return paramenters-d_paramenters*(current_itter-last_success)*10**(-O_err); #undo changes since last improvement and return new paramenters # this halts the programm
				
			else: # else try finer adjustments and rest counters
				current_itter=0
				last_success=0
				O_err=O_err+1
				return f;
				
	else: # if it didn't fail that often
		current_itter=current_itter+1
		
		if current_paramenter>=len(paramenters): # and all components are set
		
			paramenters=paramenters+d_paramenters*10**(-O_err)
			new_eval=bew(paramenters) # check the new eval 
			paramenters=paramenters+d_paramenters  # change paramenters accordingly
			last_success=last_success+isbetter(new_eval,current_eval) # increases if success, stays the same else
			current_eval=new_eval
			return f; # and repeat
			
		else: # not all components are set
			
			paramenters[current_paramenter]=paramenters[current_paramenter]+d_paramenters[current_paramenter]*10**(-O_err)*sign
			#change a component
			new_eval=bew(paramenters) # check what altering the component does
			last_success=last_success+isbetter(new_eval,current_eval) # increases if success, stays the same else
			current_eval=new_eval # and repeat
			return f;
			
			
			
			
			
			
		
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			
					stage=-1
					return optimize(paramenters,d_paramenters, current_itter, current_eval, og_eval, current_parameter_itter, stage, isbetter,current_paramenter,last_success);

