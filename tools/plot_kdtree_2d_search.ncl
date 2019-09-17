if (.not. isdefined("tag")) then
  print("Argument tag is needed!")
  exit
end if

if (.not. isdefined("ngb")) then
  print("Argument ngb is needed!")
  exit
end if
ngb = ngb - 1

wks = gsn_open_wks("pdf", "kdtree_2d")
  
node_file_path = "node_placing.node." + sprinti("%d", tag) + ".txt"
node := asciiread(node_file_path, (/numAsciiRow(node_file_path),2/), "float")

res = True
res@gsnDraw = False
res@gsnFrame = False
res@trXMinF = 0.0
res@trXMaxF = 1.0
res@trYMinF = 0.0
res@trYMaxF = 1.0

res_node = True
res_node@gsMarkerColor = "blue"
res_node@gsMarkerSizeF = 10.0

res_ngb = True
res_ngb@gsMarkerColor = "red"
res_ngb@gsMarkerSizeF = 10.0
res_ngb@gsMarkerIndex = 4

base = gsn_blank_plot(wks, res)

plot_node = gsn_add_polymarker(wks, base, node(:,0), node(:,1), res_node)

plot_ngb = gsn_add_polymarker(wks, base, node(ngb,0), node(ngb,1), res_ngb)

res_query = True
res_query@gsMarkerSizeF = 10.0
res_query@gsMarkerIndex = 12

plot_query = gsn_add_polymarker(wks, base, (/0.5/), (/0.5/), res_query)

draw(base)
frame(wks)
