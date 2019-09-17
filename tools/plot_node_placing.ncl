wks = gsn_open_wks("pdf", "node_placing")
  
do tag = 1, 380
  pdp_file_path = "node_placing.pdp." + sprinti("%d", tag) + ".txt"
  pdp := asciiread(pdp_file_path, (/numAsciiRow(pdp_file_path),2/), "float")
  
  node_file_path = "node_placing.node." + sprinti("%d", tag) + ".txt"
  node := asciiread(node_file_path, (/numAsciiRow(node_file_path),2/), "float")
  
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@trXMinF = 0.0
  res@trXMaxF = 1.0
  res@trYMinF = 0.0
  res@trYMaxF = 1.0
  
  res_pdp = True
  res_pdp@gsMarkerColor = "red"
  res_pdp@gsMarkerSizeF = 10.0
  
  res_node = True
  res_node@gsMarkerColor = "blue"
  res_node@gsMarkerSizeF = 10.0
  
  base = gsn_blank_plot(wks, res)
  
  plot_pdp = gsn_add_polymarker(wks, base, pdp(:,0), pdp(:,1), res_pdp)
  
  plot_node = gsn_add_polymarker(wks, base, node(:,0), node(:,1), res_node)
  
  draw(base)
  frame(wks)
end do
