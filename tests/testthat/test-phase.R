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

describe("phase helpers", {
  it("behave deterministically", {
    grid <- default_magic_grid(3L, 1L, step = 1L)
    expect_identical(grid, 3:1)

    result <- frsr_phase(
      phases = 4L,
      exponents = -1L:1L,
      per_cell = 3L,
      magics = grid,
      seed = 123
    )
    expect_invisible(frsr_phase_heatmap(result$heat))
  })

  it("tolerates zero or negative step sizes in magic grids", {
    zero_step <- default_magic_grid(1L, 3L, step = 0L)
    expect_identical(zero_step, 1:3)

    negative_step <- default_magic_grid(5L, -1L, step = -2L)
    expect_identical(negative_step, c(5L, 3L, 1L, -1L))
  })

  it("returns the plotted heat matrix invisibly", {
    heat <- matrix(
      c(-0.1, 0.2, 0.3, -0.4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("-2", "-1"), NULL)
    )

    returned <- frsr_phase_heatmap(heat, xlab = "Phase", ylab = "Exponent")

    expect_identical(returned, heat)
    expect_identical(rownames(returned), c("-2", "-1"))
  })
})
