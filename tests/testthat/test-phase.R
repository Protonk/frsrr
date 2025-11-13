describe("frsr_phase", {
  it("returns expected structure", {
    set.seed(11)
    result <- frsr_phase(
      phases = 8L,
      exponents = -2L:2L,
      per_cell = 4L,
      magics = c(1597413411L, 1597200000L),
      q = 0.9
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
    set.seed(99)
    good <- frsr_phase(
      phases = 4L,
      exponents = -3L:3L,
      per_cell = 6L,
      magics = 0x5f3759df
    )
    set.seed(99)
    bad <- frsr_phase(
      phases = 4L,
      exponents = -3L:3L,
      per_cell = 6L,
      magics = 0x5f100000
    )

    expect_lt(good$J, bad$J)
  })

  it("is reproducible with a fixed seed", {
    args <- list(
      phases = 6L,
      exponents = -2L:2L,
      per_cell = 4L,
      magics = c(1597413411L, 1597200000L),
      q = 0.8
    )
    set.seed(2024)
    first <- do.call(frsr_phase, args)
    set.seed(2024)
    second <- do.call(frsr_phase, args)
    expect_identical(first, second)
  })

  it("validates arguments before invoking C++ helpers", {
    expect_error(
      frsr_phase(phases = 0L),
      "`phases` must be a positive integer",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(exponents = integer()),
      "`exponents` must contain at least one value",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(exponents = 200L),
      "`exponents` must stay within [-126, 127]",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(per_cell = 0L),
      "`per_cell` must be a positive integer",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(magics = integer()),
      "`magics` must supply at least one integer constant",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(q = 0),
      "`q` must satisfy 0 < q <= 1",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(q = 2),
      "`q` must satisfy 0 < q <= 1",
      fixed = TRUE
    )
    expect_error(
      frsr_phase(NRmax = -1L),
      "`NRmax` must be a non-negative integer",
      fixed = TRUE
    )
  })
})
