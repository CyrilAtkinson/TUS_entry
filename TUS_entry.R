TUS_entry = function(target, scalp, t1=NULL, exclu=NULL, output=NULL, visual_confirm=F, resolution=0.33, voxel_size=0.1, minimal_distance=4.24, maximal_distance=7.46, transducer_size=6.4, kplan_offset=1.082) {
  
  START = Sys.time()
  
  ######################
  ### Load functions ###
  ######################
  
  require(oro.nifti)
  require(cowplot)
  require(lattice)
  require(cluster)
  
  nii_to_df = function(arr) {
    values = as.vector(arr)
    indices = expand.grid(x=1:dim(arr)[1], y=1:dim(arr)[2], z=1:dim(arr)[3])
    df = data.frame(indices, val=values)
    return(df)
  }
  
  disk_3d = function(point1, point2, radius, num_points) {
    vector_between_points = point2 - point1
    theta = seq(0, 2 * pi, length.out=num_points)
    phi = seq(0, pi, length.out=num_points)
    grid = expand.grid(theta=theta, phi=phi)
    disk_coordinates = data.frame(x=round(point1[1] + radius * cos(grid$theta) * sin(grid$phi)),
                                  y=round(point1[2] + radius * sin(grid$theta) * sin(grid$phi)),
                                  z=round(point1[3] + radius * cos(grid$phi)))
    disk_coordinates = disk_coordinates[!duplicated(disk_coordinates),]
    d = -sum(vector_between_points * point1)
    disk_coordinates$d = vector_between_points[1]*disk_coordinates$x + vector_between_points[2]*disk_coordinates$y + vector_between_points[3]*disk_coordinates$z + d
    disk_coordinates = disk_coordinates[order(abs(disk_coordinates$d)),][1:num_points,1:3]
    return(disk_coordinates)
  }
  
  sphere_3d = function(point1, radius, num_points) {
    theta = seq(0, 2 * pi, length.out=num_points)
    phi = seq(0, pi, length.out=num_points)
    grid = expand.grid(theta=theta, phi=phi)
    disk_coordinates = data.frame(x=round(point1[1] + radius * cos(grid$theta) * sin(grid$phi)),
                                  y=round(point1[2] + radius * sin(grid$theta) * sin(grid$phi)),
                                  z=round(point1[3] + radius * cos(grid$phi)))
    disk_coordinates = disk_coordinates[!duplicated(disk_coordinates),]
    return(disk_coordinates)
  }
  
  theme.novpadding = list(layout.heights = list(top.padding=0, main.key.padding=0, key.axis.padding=0, axis.xlab.padding=0, xlab.key.padding=0, key.sub.padding=0, bottom.padding=0), axis.line=list(col=0), clip=list(panel="off"), layout.widths=list(left.padding=0, key.ylab.padding=0, ylab.axis.padding=0, axis.key.padding=0, right.padding=0))
  
  ###############################
  ###        Load data        ###
  ### Verify equal dimensions ###
  ###############################
  
  empty=scalp ; empty@datatype=64 ; empty@bitpix=64 ; empty@.Data = ifelse(empty@.Data !=0, 0, 0) ; empty_data=empty@.Data
  
  if(is.null(exclu) == T){exclu=empty}
  if(is.null(t1) == T){t1=scalp}
  
  if(sd(c(dim(t1)[1], dim(target)[1], dim(scalp)[1], dim(exclu)[1])) !=0){ stop("Images provided are of different sizes") }
  if(sd(c(dim(t1)[2], dim(target)[2], dim(scalp)[2], dim(exclu)[2])) !=0){ stop("Images provided are of different sizes") }
  if(sd(c(dim(t1)[3], dim(target)[3], dim(scalp)[3], dim(exclu)[3])) !=0){ stop("Images provided are of different sizes") }
  
  ##########################
  ### Found ROI Centroid ###
  ##########################
  
  temp = nii_to_df(target@.Data)
  temp = subset(temp, val != 0)
  if(nrow(temp) > 1){
    if(nrow(temp) > 50000){temp = temp[-sample(nrow(temp) - 50000),]}
    center_target = pam(temp[1:3], 1)$medoids
  }
  if(nrow(temp) == 1){ center_target = temp[1:3] }
  rm("temp")
  
  ###################################
  ###           Step 1            ###
  ### Found possible entry points ###
  ###      Based on distance      ###
  ###################################
  
  for(i in 1:dim(scalp)[1]){ if (sum(scalp@.Data[i, 1:dim(scalp)[2], 1:dim(scalp)[3]], na.rm=T) != 0){xmin=i-1 ; break}} ; if(xmin==0){xmin=1}
  for(i in dim(scalp)[1]:1){ if (sum(scalp@.Data[i, 1:dim(scalp)[2], 1:dim(scalp)[3]], na.rm=T) != 0){xmax=i+1 ; break}} ; if(xmax>dim(scalp)[1]){xmax=dim(scalp)[1]}
  for(i in 1:dim(scalp)[2]){ if (sum(scalp@.Data[1:dim(scalp)[1], i, 1:dim(scalp)[3]], na.rm=T) != 0){ymin=i-1 ; break}} ; if(ymin==0){ymin=1}
  for(i in dim(scalp)[2]:1){ if (sum(scalp@.Data[1:dim(scalp)[1], i, 1:dim(scalp)[3]], na.rm=T) != 0){ymax=i+1 ; break}} ; if(ymax>dim(scalp)[2]){ymax=dim(scalp)[2]}
  for(i in 1:dim(scalp)[3]){ if (sum(scalp@.Data[1:dim(scalp)[1], 1:dim(scalp)[2], i], na.rm=T) != 0){zmin=i-1 ; break}} ; if(zmin==0){zmin=1}
  for(i in dim(scalp)[3]:1){ if (sum(scalp@.Data[1:dim(scalp)[1], 1:dim(scalp)[2], i], na.rm=T) != 0){zmax=i+1 ; break}} ; if(zmax>dim(scalp)[3]){zmax=dim(scalp)[3]}
  
  delta = round((maximal_distance / voxel_size) +1)
  perimetre = round(2 * pi * delta)
  box = sphere_3d(point1=center_target, radius=delta, num_points=perimetre)
  box$x = ifelse(box$x<xmin, xmin, ifelse(box$x>xmax, xmax, box$x))
  box$y = ifelse(box$y<ymin, ymin, ifelse(box$y>ymax, ymax, box$y))
  box$z = ifelse(box$z<zmin, zmin, ifelse(box$z>zmax, zmax, box$z))
  box = box[!duplicated(box),]
  box = box[scalp[cbind(c(box$x), c(box$y), c(box$z))] == 0,]
  
  vect = max(abs(cbind(c(center_target[1]-box$x), c(center_target[2]-box$y), c(center_target[3]-box$z))))
  
  rm("i","xmin","xmax","ymin","ymax","zmin","zmax","delta","perimetre")
  
  Final_1_distance = empty ; Final_data_1_distance = empty_data
  
  NOW = Sys.time()
  
  for(i in 1:nrow(box)){
    if(difftime(Sys.time(), NOW, units="secs") > 1){
      difference_time = difftime(Sys.time(), START, units="secs")
      time_min = trunc(difference_time/60)
      time_sec = trunc(difference_time - (time_min*60))
      NOW = Sys.time()
      cat(paste("\014", output,
                "\nStep 1/4: Distance selection       - ", round((i/nrow(box))*100,1), "%",
                "\nStep 2/4: Decrease resolution      - wait",
                "\nStep 3/4: Angle calculation        - wait",
                "\nStep 4/4: Calculate exact distance - wait",
                "\nTotal duration [min]: ", time_min, ".", time_sec, sep=""))
    }
    vect_mat = cbind(c(round(seq(box$x[i], center_target[1], length=vect))), c(round(seq(box$y[i], center_target[2], length=vect))), c(round(seq(box$z[i], center_target[3], length=vect))))
    if(sum(exclu[vect_mat]) != 0){next}
    vect_mat = vect_mat[scalp[vect_mat] == 0,]
    if(isTRUE(nrow(vect_mat) > 1)){vect_mat = vect_mat[nrow(vect_mat),]}
    if(length(vect_mat) > 0){
      Dist = dist(matrix(c(center_target, vect_mat), nrow=2, byrow=T))*voxel_size
      if(Dist <= maximal_distance){
        Final_data_1_distance[vect_mat[1],vect_mat[2],vect_mat[3]] = Dist
      }
    }
  }
  
  Final_1_distance@.Data = Final_data_1_distance
  
  rm("box","Dist","i","vect","vect_mat","difference_time","time_min","time_sec","NOW")
  
  ###########################
  ###       Step 2        ###
  ### Decrease resolution ###
  ###    to save time     ###
  ###########################
  
  data = nii_to_df(Final_1_distance) ; data = subset(data, val != 0) ; minimal_depth = min(data$val, na.rm=T)
  data_new = data ; data_new$Done = 0
  
  if(nrow(data) == 0){ stop("Target not reachable (too deep)") }
  
  Final_2_resolution = empty ; Final_data_2_resolution = empty_data
  
  NOW = Sys.time()
  
  if(resolution > voxel_size){
    for(i in 1:nrow(data)){
      if(i >= nrow(data_new)){break}
      if(difftime(Sys.time(), NOW, units="secs") > 1){
        difference_time = difftime(Sys.time(), START, units="secs")
        time_min = trunc(difference_time/60)
        time_sec = trunc(difference_time - (time_min*60))
        NOW = Sys.time()
        cat(paste("\014", output,
                  "\nStep 1/4: Distance selection       - done",
                  "\nStep 2/4: Decrease resolution      - in progress",
                  "\nStep 3/4: Angle calculation        - wait",
                  "\nStep 4/4: Calculate exact distance - wait",
                  "\nTotal duration [min]: ", time_min, ".", time_sec, sep=""))
      }
      
      data_new$Done[i] = "OK"
      data_new$Done = ifelse(data_new$Done == "OK", "OK", ifelse(
        (data_new$x >= data_new$x[i]-resolution/voxel_size & data_new$x <= data_new$x[i]+resolution/voxel_size &
           data_new$y >= data_new$y[i]-resolution/voxel_size & data_new$y <= data_new$y[i]+resolution/voxel_size &
           data_new$z >= data_new$z[i]-resolution/voxel_size & data_new$z <= data_new$z[i]+resolution/voxel_size), "Close", "Far"))
      data_new_temp = data_new[data_new$Done == "Close",]
      if(nrow(data_new_temp) == 0){next}
      data_new_temp$x2 = data_new$x[i] ; data_new_temp$y2 = data_new$y[i] ; data_new_temp$z2 = data_new$z[i]
      data_new_temp$resolution_distance = apply(data_new_temp[,c(1:3,6:8)], 1, function(x) dist(matrix(x, nrow=2, byrow=T)))*voxel_size
      data_new_temp = data_new_temp[data_new_temp$resolution_distance >= resolution,]
      data_new = rbind(data_new[data_new$Done == "OK",], data_new_temp[,c(1:5)], data_new[data_new$Done == "Far",])
    }
    Final_data_2_resolution[cbind(c(data_new$x), c(data_new$y), c(data_new$z))] = data_new$val
    
  }
  
  if(resolution <= voxel_size){Final_data_2_resolution = Final_data_1_distance}
  
  Final_2_resolution@.Data = Final_data_2_resolution
  
  rm("data","data_new","data_new_temp","i","difference_time","time_min","time_sec","NOW")
  
  #########################
  ###      Step 3       ###
  ### Angle calculation ###
  #########################
  
  data = nii_to_df(Final_2_resolution) ; data = subset(data, val != 0)
  
  Final_3_angle = empty ; Final_data_3_angle = empty_data
  
  delta = (transducer_size / 2 ) / voxel_size ; perimetre = round(2 * pi * delta)
  
  NOW = Sys.time()
  
  for (i in 1:nrow(data)){
    exclusion = "n"
    if(difftime(Sys.time(), NOW, units="secs") > 1){
      difference_time = difftime(Sys.time(), START, units="secs")
      time_min = trunc(difference_time/60)
      time_sec = trunc(difference_time - (time_min*60))
      NOW = Sys.time()
      cat(paste("\014", output,
                "\nStep 1/4: Distance selection       - done",
                "\nStep 2/4: Decrease resolution      - done",
                "\nStep 3/4: Angle calculation        - ", round((i/nrow(data))*100,1), "%",
                "\nStep 4/4: Calculate exact distance - wait",
                "\nTotal duration [min]: ", time_min, ".", time_sec, sep=""))
    }
    
    disk = disk_3d(point1=c(data[i,1], data[i,2], data[i,3]), point2=center_target, radius=delta, num_points = perimetre)
    disk_in = disk[disk$x>0 & disk$y>0 & disk$z>0 & disk$x<=dim(scalp)[1] & disk$y<=dim(scalp)[2] & disk$z<=dim(scalp)[3],]
    disk_in$loc = scalp[cbind(c(disk_in$x), c(disk_in$y), c(disk_in$z))]
    disk_out = rbind(disk[disk$x<=0 | disk$y<=0 | disk$z<=0 | disk$x>dim(scalp)[1] | disk$y>dim(scalp)[2] | disk$z>dim(scalp)[3],], disk_in[disk_in$loc == 0,1:3])
    disk_in = disk_in[disk_in$loc == 1,]
    
    if(nrow(disk_out) != 0){
      
      vect_out = max(abs(cbind(disk_out$x-center_target[1], disk_out$y-center_target[2], disk_out$z-center_target[3])))
      disk_out$dist = 0
      
      for(ii in 1:nrow(disk_out)){
        vect_mat = cbind(c(round(seq(center_target[1], disk_out$x[ii], length=vect_out))), c(round(seq(center_target[2], disk_out$y[ii], length=vect_out))), c(round(seq(center_target[3], disk_out$z[ii], length=vect_out))))
        vect_mat = vect_mat[vect_mat[,1]>0 & vect_mat[,2]>0 & vect_mat[,3]>0 & vect_mat[,1]<=dim(scalp)[1] & vect_mat[,2]<=dim(scalp)[2] & vect_mat[,3]<=dim(scalp)[3],]
        if(sum(exclu[cbind(c(vect_mat[,1]), c(vect_mat[,2]), c(vect_mat[,3]))]) != 0){exclusion = "y" ; break}
        vect_mat = vect_mat[scalp[cbind(c(vect_mat[,1]), c(vect_mat[,2]), c(vect_mat[,3]))] == 0,]
        if(length(vect_mat) > 3){vect_mat = vect_mat[1,]}
        if(length(vect_mat) > 0){ disk_out$dist[ii] = dist(matrix(c(disk_out$x[ii], disk_out$y[ii], disk_out$z[ii], vect_mat), nrow=2, byrow=T))*voxel_size }
      }
    }
    
    if(exclusion == "y"){next}
    
    if(nrow(disk_in) != 0){
      
      vect_in = cbind(disk_in$x-center_target[1], disk_in$y-center_target[2], disk_in$z-center_target[3])
      vect_in = disk_in + (vect_in * 2)
      vect_in_max = max(abs(vect_in))
      disk_in$dist = 0
      
      for(ii in 1:nrow(disk_in)){
        vect_mat = cbind(c(round(seq(center_target[1], vect_in$x[ii], length=vect_in_max))), c(round(seq(center_target[2], vect_in$y[ii], length=vect_in_max))), c(round(seq(center_target[3], vect_in$z[ii], length=vect_in_max))))
        vect_mat = vect_mat[vect_mat[,1]>0 & vect_mat[,2]>0 & vect_mat[,3]>0 & vect_mat[,1]<=dim(scalp)[1] & vect_mat[,2]<=dim(scalp)[2] & vect_mat[,3]<=dim(scalp)[3],]
        if(sum(exclu[cbind(c(vect_mat[,1]), c(vect_mat[,2]), c(vect_mat[,3]))]) != 0){exclusion = "y" ; break}
        vect_mat = vect_mat[scalp[cbind(c(vect_mat[,1]), c(vect_mat[,2]), c(vect_mat[,3]))] == 0,]
        if(length(vect_mat) > 3){vect_mat = vect_mat[1,]}
        if(length(vect_mat) > 0){ disk_in$dist[ii] = dist(matrix(c(disk_in$x[ii], disk_in$y[ii], disk_in$z[ii], vect_mat), nrow=2, byrow=T))*voxel_size }
      }
    }
    
    if(exclusion == "y"){next}
    
    if(nrow(disk_out) != 0 & nrow(disk_in) != 0){Final_data_3_angle[data[i,1], data[i,2], data[i,3]] = sd(c(-disk_in$dist, disk_out$dist), na.rm=T)}
    if(nrow(disk_out) == 0 & nrow(disk_in) != 0){Final_data_3_angle[data[i,1], data[i,2], data[i,3]] = sd(c(-disk_in$dist), na.rm=T)}
    if(nrow(disk_out) != 0 & nrow(disk_in) == 0){Final_data_3_angle[data[i,1], data[i,2], data[i,3]] = sd(c(disk_out$dist), na.rm=T)}
    
  }
  
  Final_3_angle@.Data = Final_data_3_angle
  
  rm("data","disk","disk_in","disk_out","vect_mat","vect_in","vect_out","exclusion","vect_in_max","delta","difference_time","i","ii","perimetre","time_min","time_sec","NOW")
  
  #############################
  ###        Step 4         ###
  ### Verify exact distance ###
  ###   Calculate volume    ###
  #############################
  
  data = nii_to_df(Final_3_angle) ; data = subset(data, val != 0) ; data = data[order(data$val),]
  
  if(nrow(data) == 0){ stop("Target not reachable (exclusion area on the way)") }
  
  delta = (transducer_size / 2 ) / voxel_size ; perimetre = round(2 * pi * delta)
  
  Final_T1 = t1
  happy=1
  
  for (i in 1:nrow(data)){
    difference_time = difftime(Sys.time(), START, units="secs")
    time_min = trunc(difference_time/60)
    time_sec = trunc(difference_time - (time_min*60))
    cat(paste("\014", output,
              "\nStep 1/4: Distance selection       - done",
              "\nStep 2/4: Decrease resolution      - done",
              "\nStep 3/4: Angle calculation        - done",
              "\nStep 4/4: Calculate exact distance - ", i, "/", nrow(data),
              "\nTotal duration [min]: ", time_min, ".", time_sec, sep=""))
    
    # found the right position
    
    vect = data.frame(x=data$x[i]-center_target[1], y=data$y[i]-center_target[2], z=data$z[i]-center_target[3])
    vect_diff = vect/max(abs(vect))
    
    out = "no" ; iter = 0
    while(out == "no"){
      x1 = round(data[i,1]+(vect_diff[1]*iter)[[1]])
      y1 = round(data[i,2]+(vect_diff[2]*iter)[[1]])
      z1 = round(data[i,3]+(vect_diff[3]*iter)[[1]])
      disk = disk_3d(point1 = c(x1, y1, z1), point2 = center_target, radius = delta, num_points = perimetre)[,1:3]
      disk_select = disk[disk$x>0 & disk$y>0 & disk$z>0 & disk$x<=dim(scalp)[1] & disk$y<=dim(scalp)[2] & disk$z<=dim(scalp)[3],]
      if(sum(scalp[cbind(c(disk_select$x), c(disk_select$y), c(disk_select$z))]) == 0){break}
      if(nrow(disk_select) == 0){break}
      iter = iter+1
    }
    
    # calculate the right distance and take another point if out of the range
    
    final_distance = round(dist(matrix(c(x1,y1,z1, center_target), nrow=2, byrow=T))*voxel_size,3)
    if((final_distance<minimal_distance | final_distance>maximal_distance) & i == nrow(data)){ stop("Target not reachable (nothing survived to criteria)") }
    if(final_distance<minimal_distance | final_distance>maximal_distance){next}
    
    ## fill the disk
    
    vect = max(abs(cbind(c(x1-disk$x), c(y1-disk$y), c(z1-disk$z))))
    for(ii in 1:nrow(disk)){
      vect_mat = cbind(c(round(seq(disk$x[ii], x1, length=vect))), c(round(seq(disk$y[ii], y1, length=vect))), c(round(seq(disk$z[ii], z1, length=vect))))
      if(ii == 1){disk_fill = vect_mat}
      if(ii != 1){disk_fill = rbind(disk_fill, vect_mat)}
    }
    disk_fill = data.frame(disk_fill[!duplicated(disk_fill),])
    colnames(disk_fill) = c("x","y","z")
    
    ## Generate ROI data
    
    vect = max(abs(cbind(c(center_target[1]-x1), c(center_target[2]-y1), c(center_target[3]-z1))))
    vect = data.frame(x=round(seq(x1, center_target[1], length=vect)), y=round(seq(y1, center_target[2], length=vect)), z=round(seq(z1, center_target[3], length=vect)))
    
    Marks = rbind(
      data.frame(expand.grid(x=seq(center_target[1]-1,center_target[1]+1), y=seq(center_target[2]-1,center_target[2]+1), z=seq(center_target[3]-1,center_target[3]+1)), cat="target"),
      data.frame(expand.grid(x=seq(data$x[i]-1,data$x[i]+1), y=seq(data$y[i]-1,data$y[i]+1), z=seq(data$z[i]-1,data$z[i]+1)), cat="entry"),
      data.frame(expand.grid(x=seq(x1-1,x1+1), y=seq(y1-1,y1+1), z=seq(z1-1,z1+1)), cat="transducer"),
      data.frame(vect, cat="vector"),
      data.frame(disk_fill, cat="disk")
    )
    
    Marks = subset(Marks, x>0 & y>0 & z>0 & x<=dim(scalp)[1] & y<=dim(scalp)[2] & z<=dim(scalp)[3])
    
    ## Visual approvement if required
    
    if(isTRUE(visual_confirm)){
      Final_data_T1 = t1@.Data
      if(max(t1, na.rm=T) <= 1){Final_data_T1[cbind(c(Marks$x), c(Marks$y), c(Marks$z))] = 2}
      if(max(t1, na.rm=T) > 1){Final_data_T1[cbind(c(Marks$x), c(Marks$y), c(Marks$z))] = max(Final_data_T1, na.rm=T)}
      p1 = levelplot(Final_data_T1[data$x[i],,], colorkey=F, xlab="", ylab="", scales=list(x=list(draw=F), y=list(draw=F)),par.settings=theme.novpadding)
      p2 = levelplot(Final_data_T1[,data$y[i],], colorkey=F, xlab="", ylab="", scales=list(x=list(draw=F), y=list(draw=F)),par.settings=theme.novpadding)
      p3 = levelplot(Final_data_T1[,,data$z[i]], colorkey=F, xlab="", ylab="", scales=list(x=list(draw=F), y=list(draw=F)),par.settings=theme.novpadding)
      plot(plot_grid(p1, p2, p3, labels=paste(i, "/", nrow(data), "- Angle =",round(data$val[i],1))))
      Sys.sleep(1)
      CHOICE = askYesNo("Are you happy with this transducer position?") 
      if (CHOICE == "FALSE"){ graphics.off() ; happy = 0 }
      if (CHOICE == "TRUE"){ happy = 1 }
    }
    
    if(happy == 1){break}
    
  }
  
  ## Calculate volume between transducer and head
  
  
  vect = max(abs(cbind(c(center_target[1]-disk$x), c(center_target[2]-disk$y), c(z=center_target[3]-disk$z))))
  for(ii in 1:nrow(disk_fill)){
    vect_mat = cbind(c(round(seq(disk_fill$x[ii], center_target[1], length=vect))), c(round(seq(disk_fill$y[ii], center_target[2], length=vect))), c(round(seq(disk_fill$z[ii], center_target[3], length=vect))))[-1,]
    vect_mat_out = vect_mat[vect_mat[,1]<=0 | vect_mat[,2]<=0 | vect_mat[,3]<=0 | vect_mat[,1]>dim(scalp)[1] | vect_mat[,2]>dim(scalp)[2] | vect_mat[,3]>dim(scalp)[3],]
    vect_mat_in = vect_mat[vect_mat[,1]>0 & vect_mat[,2]>0 & vect_mat[,3]>0 & vect_mat[,1]<=dim(scalp)[1] & vect_mat[,2]<=dim(scalp)[2] & vect_mat[,3]<=dim(scalp)[3],]
    vect_mat_out = rbind(vect_mat_out, vect_mat_in[scalp[cbind(c(vect_mat_in[,1]), c(vect_mat_in[,2]), c(vect_mat_in[,3]))] == 0,])
    if(ii == 1){volume = vect_mat_out}
    if(ii != 1){volume = rbind(volume, vect_mat_out)}
  }
  volume = data.frame(volume[!duplicated(volume),])
  colnames(volume) = c("x","y","z")
  
  ############
  ### Save ###
  ############
  
  Final_T1 = t1
  Marks_neuronav = subset(Marks, cat == "entry" | cat == "target")
  if(max(t1, na.rm=T) <= 1){Final_T1@.Data[cbind(c(Marks_neuronav$x), c(Marks_neuronav$y), c(Marks_neuronav$z))] = 2}
  if(max(t1, na.rm=T) > 1){Final_T1@.Data[cbind(c(Marks_neuronav$x), c(Marks_neuronav$y), c(Marks_neuronav$z))] = max(Final_T1, na.rm=T)}
  
  Final_3D = empty ; Final_3D@.Data = empty_data
  Final_3D@.Data[cbind(c(Marks$x), c(Marks$y), c(Marks$z))] = 1
  
  Final_ROI = empty ; Final_ROI@.Data = empty_data
  Marks_entry = subset(Marks, cat == "entry")
  Final_ROI@.Data[cbind(c(Marks_entry$x), c(Marks_entry$y), c(Marks_entry$z))] = 1
  
  Final_volume = empty ; Final_volume@.Data = empty_data
  volume_in = subset(volume, x>0 & y>0 & z>0 & x<=dim(t1)[1] & y<=dim(t1)[2] & z<=dim(t1)[3])
  Final_volume@.Data[cbind(c(volume_in$x), c(volume_in$y), c(volume_in$z))] = 1
  
  Difference_time = difftime(Sys.time(), START, units="secs")
  time_min = trunc(Difference_time/60)
  time_sec = trunc(Difference_time - (time_min*60))
  
  report = paste(
    "ID: ", output,
    "\n\nDistance: ", final_distance*10, " mm",
    "\nAngle score: ", round(data$val[i], 4),
    "\nVolume: ", nrow(volume), " voxels",
    "\n\nTarget coordinates [world]:     ", round(translateCoordinate(c(center_target), t1)[1]), " . ", round(translateCoordinate(c(center_target), t1)[2]),  " . ", round(translateCoordinate(c(center_target), t1)[3]),
    "\nTransducer coordinates [world]: ", round(translateCoordinate(x1, t1)[1]), " . ", round(translateCoordinate(y1, t1)[2]),  " . ", round(translateCoordinate(z1, t1)[3]),
    "\n\nDistance [kplan]: ", (final_distance*10) - (kplan_offset*10), " mm",
    "\nTarget coordinates [kplan]:     ", dim(Final_T1)[1]-center_target[1], " . ", dim(Final_T1)[2]-center_target[2], " . ", center_target[3]-1,
    "\nTransducer coordinates [kplan]: ", dim(Final_T1)[1]-x1, " . ", dim(Final_T1)[2]-y1, " . ", z1-1,
    "\n\nComputation time: ", time_min, ".", time_sec, " mins",
    sep="")
  
  list_return = list("report"=report,
                     "target"=target, "t1"=t1, "scalp"=scalp, "exclu"=exclu, "output"=output, "visual_confirm"=visual_confirm, "resolution"=resolution, "voxel_size"=voxel_size, "minimal_distance"=minimal_distance, "maximal_distance"=maximal_distance, "transducer_size"=transducer_size, "kplan_offset"=kplan_offset,
                     "Target_coordinates"=center_target, "Scalp_coordinates"=data[i,1:3], "Transducer_coordinates"=c(x1,y1,z1), "Distance"=final_distance, "Angle"=data$val[i], "Volume_k"=nrow(volume), "minimal_depth"=minimal_depth, "computation_time"=paste(time_min, ".", time_sec, " mins", sep=""), "MRI_1_distance"=Final_1_distance, "MRI_2_resolution"=Final_2_resolution, "MRI_3_Angle"=Final_3_angle, "MRI_Final_neuronav"=Final_T1, "MRI_Final_validation"=Final_3D, "MRI_Final_ROI"=Final_ROI, "MRI_Final_volume"=Final_volume
  )
  
  return(list_return)
  
}



