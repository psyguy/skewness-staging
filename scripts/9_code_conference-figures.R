ts.length <- 5000


d.3.nar <- make_dgm_df(
  dgm_nar(
    phi = 0.4,
    var.resid = 47,
    T = ts.length,
    Mean = 85,
    seed = 1 + 9
  ),
  dgm_nar(
    phi = 0.4,
    var.resid = 20,
    T = ts.length,
    Mean = 45,
    seed = 2 + 4
  ),
  dgm_nar(
    phi = 0.4,
    var.resid = 10,
    T = ts.length,
    Mean = 20,
    seed = 3 + 2
  )
)


profile_nar <- plot_dgm.profile(d.3.nar,
                                brewer.pal(name = "YlOrBr", n = 9)[c(5, 7, 9)],
                                "$AR(1)$")


ggsave(
  "iops-nar-0.4.pdf",
  profile_nar,
  path = "figures",
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)




d.3.binar <- make_dgm_df(
  dgm_generator(
    Model = "BinAR",
    k = 7,
    Mean = 3.5,
    phi = 0.4,
    T = ts.length,
    seed = 1 + 5
  ),
  dgm_generator(
    Model = "BinAR",
    k = 7,
    Mean = 2,
    phi = 0.4,
    T = ts.length,
    seed = 2 + 6
  ),
  dgm_generator(
    Model = "BinAR",
    k = 7,
    Mean = 0.5,
    phi = 0.4,
    T = ts.length,
    seed = 3
  )
)

profile_binar <- plot_dgm.profile(d.3.binar,
                                  brewer.pal(name = "YlGn", n = 9)[c(5, 7, 9)],
                                  "$BinAR(1)$")


ggsave(
  "iops-binar-0.4.pdf",
  profile_binar,
  path = "figures",
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)





d.3.podar <- make_dgm_df(
  dgm_generator(
    Model = "PoDAR",
    tau = 0.4,
    Mean = 40,
    T = ts.length,
    seed = 1 + 4
  ),
  dgm_generator(
    Model = "PoDAR",
    tau = 0.4,
    Mean = 10,
    T = ts.length,
    seed = 2 + 4
  ),
  dgm_generator(
    Model = "PoDAR",
    tau = 0.4,
    Mean = 1.,
    T = ts.length,
    seed = 3 + 2
  )
)


profile_podar <- plot_dgm.profile(d.3.podar,
                                  brewer.pal(name = "BuPu", n = 9)[c(5, 7, 9)],
                                  "$\\PoDAR(1)$")


ggsave(
  "iops-podar-0.4.pdf",
  profile_podar,
  path = "figures",
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)



d.3.chiar <- make_dgm_df(
  dgm_generator(
    Model = "ChiAR",
    phi = 0.4,
    nu = 25,
    T = ts.length,
    seed = 1 + 1
  ),
  dgm_generator(
    Model = "ChiAR",
    phi = 0.4,
    nu = 5,
    T = ts.length,
    seed = 2 + 5
  ),
  dgm_generator(
    Model = "ChiAR",
    phi = 0.4,
    nu = 1,
    T = ts.length,
    seed = 3 + 1
  )
)


profile_chiar <- plot_dgm.profile(d.3.chiar,
                                  brewer.pal(name = "GnBu", n = 9)[c(5, 7, 9)],
                                  "$\\chi^2AR(1)$")



ggsave(
  "iops-chiar-0.4.pdf",
  profile_chiar,
  path = "figures",
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)


