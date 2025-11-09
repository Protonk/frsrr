describe("frsr_phase", {
  it("returns expected structure", {
    result <- frsr_phase(
      phases = 8L,
      exponents = -2L:2L,
      per_cell = 4L,
      magics = c(1597413411L, 1597200000L),
      q = 0.9,
      seed = 11
    )

    expect_type(result$magic, "integer")
    expect_type(result$J, "double")
    expect_type(result$R, "double")
    expect_s3_class(result$phase_tbl, "data.frame")
    expect_true(is.matrix(result$heat))
    expect_equal(nrow(result$phase_tbl), 8)
    expect_equal(result$phase_tbl$phase_id, seq_len(8))
    expect_true(all(result$phase_tbl$n == length(-2L:2L) * 4L))
    expect_equal(nrow(result$heat), length(-2L:2L))
    expect_equal(ncol(result$heat), 8)
  })

  it("worsens objective for clearly bad magic constants", {
    good <- frsr_phase(
      phases = 4L,
      exponents = -3L:3L,
      per_cell = 6L,
      magics = 0x5f3759df,
      seed = 99
    )
    bad <- frsr_phase(
      phases = 4L,
      exponents = -3L:3L,
      per_cell = 6L,
      magics = 0x5f100000,
      seed = 99
    )

    expect_lt(good$J, bad$J)
  })
})
