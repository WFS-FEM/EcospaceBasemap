#' Input data contract for the WFS Ecospace basemaps.
#'
#' Most of data/ is gitignored -- roughly 628 MB, and two of the inputs cannot be
#' redistributed at all. A clone therefore starts with 3 of the 9 inputs present
#' and the rest have to be downloaded or requested. This file is the single place
#' that records what those inputs are, where they come from, and what shape they
#' have to be in on disk, so the driver, its error messages and the README cannot
#' drift apart.
#'
#' Three entry points:
#'   fn.data_manifest()  the contract itself, as a data.frame
#'   fn.check_inputs()   what is present, what is missing, and how to get it
#'   fn.pull_all()       run every automatic download that is still needed
#'
#' base R only -- this has to work before any spatial package is attached.


#' The input data contract.
#'
#' One row per input the driver needs. `path` is relative to `dir.data` (the
#' repo's data/ folder unless overridden in config.local.R).
#'
#' Columns:
#'   key       short identifier used in code and messages
#'   section   section of make_WFS_basemaps.R that consumes it
#'   path      location relative to dir.data ("" for the geodatabase, which is
#'             discovered by extension)
#'   type      how presence is tested: dir.shp, dir.asc, file, gdb, dir.any
#'   required  TRUE stops the run when missing; FALSE degrades with a warning
#'   how       repo (ships with the clone), auto (downloads itself), manual
#'   pull      name of the function that fetches it, NA for repo and manual
#'   pull.dir  directory, relative to dir.data, that the puller is handed. Held
#'             explicitly rather than derived from `path`, because the pullers
#'             differ in what they create: fn.pull_dbseabed() makes its own
#'             subdirectories inside what it is given, while fn.pull_seagrass()
#'             is handed the parent of the directory it creates
#'   url       direct download, or the page to obtain it from
#'   terms     redistribution note
#'   expects   the shape the pipeline requires on disk
#'
#' @return data.frame, one row per input.
fn.data_manifest <- function() {

  r <- function(key, section, path, type, required, how, pull, pull.dir,
                url, terms, expects)
    data.frame(key = key, section = section, path = path, type = type,
               required = required, how = how, pull = pull, pull.dir = pull.dir,
               url = url, terms = terms, expects = expects,
               stringsAsFactors = FALSE)

  rbind(
    r("regions", "5", "regions", "dir.any", TRUE, "repo", NA_character_, NA_character_,
      "ships with the repository",
      "in the repository",
      paste("age0_survey_regions.shp (+ .dbf/.prj/.shx), env3LABS_93to24.csv, and the",
            "HAND-EDITED age0_survey_regions_5min_mod.asc, which is an INPUT and exists",
            "nowhere else -- never overwrite it with generated output")),

    r("management_areas", "3", "management_areas", "dir.any", TRUE, "repo", NA_character_, NA_character_,
      "ships with the repository",
      "in the repository (originally Gulf Council / NOAA SERO)",
      "8 zipped shapefiles, one per management area"),

    r("ports", "4", "ports", "dir.any", TRUE, "repo", NA_character_, NA_character_,
      "ships with the repository",
      "in the repository (originally FWC ReportCreator; NOAA MRIP)",
      paste("ReportCreatorResults-County.csv, from a specific 1985-2025 all-species",
            "ReportCreator query, and MRIP WFS gag and red grouper dtrips by county.csv")),

    r("dbseabed", "2.2", "dbseabed", "dir.asc", TRUE, "auto", "fn.pull_dbseabed", "dbseabed",
      "https://csdms.colorado.edu/wiki/DBSEABED#Data_for_Modellers",
      "public",
      "four subdirectories Gmf_GVL/ Gmf_MUD/ Gmf_RCK/ Gmf_SND/, each holding .asc grids"),

    r("seagrass_gulfwide", "2.1", "seagrass/GulfwideSAV", "dir.shp", TRUE, "auto",
      "fn.pull_seagrass", "seagrass",
      "https://www.ncei.noaa.gov/waf/data-atlas-waf/biotic/documents/GulfwideSAV.zip",
      "public (NOAA NCEI Gulf Data Atlas)",
      "a directory holding one .shp; roughly 90k polygons, whole-Gulf extent"),

    r("seagrass_fwc", "2.1", "seagrass/Seagrass_Statewide", "dir.shp", FALSE, "auto",
      "fn.pull_seagrass_fwc", "seagrass",
      "https://geodata.myfwc.com/datasets/myfwc::seagrass-habitat-in-florida",
      "public (FWC open data)",
      paste("a directory named exactly Seagrass_Statewide holding one .shp -- the name",
            "is used in output filenames and plot titles. Roughly 230 MB. Optional:",
            "without it section 2.1 falls back to GulfwideSAV alone")),

    r("artificial_reefs_fwc", "2.4", "artificial_reefs/reeflocations.csv", "file", FALSE,
      "auto", "fn.pull_reeflocations", "artificial_reefs",
      "https://geodata.myfwc.com/datasets/artificial-reefs-in-florida",
      "public (FWC open data)",
      paste("CSV with a leading row-name column plus Description (free text) and Relief;",
            "read with read.csv(row.names = 1, fileEncoding = 'windows-1252').",
            "Optional: without it section 2.4 runs unweighted with no relief classes")),

    r("gfisher_gdb", "2.3", "", "gdb", TRUE, "manual", NA_character_, NA_character_,
      "supplied by FWRI on request -- not publicly downloadable",
      paste("FWRI-supplied; check redistribution terms before sharing. 274 MB, so it",
            "cannot be committed"),
      paste("a *.gdb directory placed directly in dir.data, e.g.",
            "GFISHER_EAST_Universe_2026.gdb. Discovered by extension, or named",
            "explicitly with file.gdb in config.local.R")),

    r("artificial_reefs_structures", "2.4",
      "artificial_reefs/dataS2_artificial_reef_structures_REDACTED.csv",
      "file", TRUE, "manual", NA_character_, NA_character_,
      "published artificial reef structure database (supplementary data table S2)",
      "marked REDACTED by its source; check terms before redistributing",
      paste("CSV with a leading row-name column and the columns state, description,",
            "lat_dd, long_dd, area_m2; the FL_GOM rows are the ones used"))
  )
}


#' Reference item counts, for detecting that an input has changed underneath you.
#'
#' Measured on the copies that produced the verified 5 arc-min basemaps -- the run
#' whose habitat figure came out identical to docs/img/02-5-habitat-sum1.png. Units
#' follow fn.input_unit(): features for shapefiles, LINES for CSVs (not records --
#' see fn.input_fingerprint()), files otherwise.
#'
#' These are a tripwire, not a specification. FWC revises the seagrass compilation
#' and grows the artificial reef table, so a difference here is often legitimate;
#' it just should never pass unnoticed. The geodatabase entry earns its place for a
#' different reason: 79 files is how you learn that a 274 MB copy truncated.
#'
#' Update a value only alongside a run that shows the new data still produces a
#' sensible basemap, and say so in the commit.
#'
#' @return Named integer vector, one element per manifest key. NA disables the
#'   check for that input.
fn.reference_counts <- function()
  c(regions                     = 6L,       # shapefile + sidecars, .asc, .csv
    management_areas            = 8L,       # one zip per area
    ports                       = 2L,       # FWC ReportCreator + NOAA MRIP
    dbseabed                    = 4L,       # Gmf_GVL/MUD/RCK/SND
    seagrass_gulfwide           = 30335L,
    seagrass_fwc                = 86173L,   # FWC compilation as of Sept 2026
    artificial_reefs_fwc        = 4550L,    # 4,548 records; 2 hold embedded newlines
    gfisher_gdb                 = 79L,
    artificial_reefs_structures = 14428L)


#' Test whether one manifest row is satisfied on disk.
#'
#' @param row One row of fn.data_manifest().
#' @param dir.data Data root.
#' @param file.gdb Explicit geodatabase path, or NULL to discover one.
#' @return list(present, found); found is the resolved path, or NA.
fn.input_found <- function(row, dir.data, file.gdb = NULL) {

  p <- if (nzchar(row$path)) file.path(dir.data, row$path) else dir.data

  switch(row$type,

    "file" = list(present = file.exists(p),
                  found   = if (file.exists(p)) p else NA_character_),

    "dir.shp" = {
      f <- if (dir.exists(p)) list.files(p, pattern = "\\.shp$", full.names = TRUE) else character(0)
      list(present = length(f) > 0, found = if (length(f) > 0) f[1] else NA_character_)
    },

    "dir.asc" = {
      f <- if (dir.exists(p)) list.files(p, pattern = "\\.asc$", recursive = TRUE) else character(0)
      list(present = length(f) >= 4, found = if (length(f) >= 4) p else NA_character_)
    },

    "dir.any" = {
      f <- if (dir.exists(p)) list.files(p) else character(0)
      list(present = length(f) > 0, found = if (length(f) > 0) p else NA_character_)
    },

    "gdb" = {
      g <- if (!is.null(file.gdb) && nzchar(file.gdb)) file.gdb
           else list.files(dir.data, pattern = "\\.gdb$", full.names = TRUE)
      g <- g[dir.exists(g) | file.exists(g)]
      list(present = length(g) > 0, found = if (length(g) > 0) g[1] else NA_character_)
    },

    stop("Unknown manifest type: ", row$type)
  )
}


#' Human-readable unit for an input's item count.
#'
#' Derived from `type` rather than stored per row, since it never varies within a
#' type. "lines" is deliberate for CSVs: see fn.input_fingerprint().
#'
#' @param type A manifest `type` value.
#' @return Single string.
fn.input_unit <- function(type)
  switch(type, dir.shp = "features", file = "lines", dir.asc = "grids", "files")


#' Fingerprint one input as it currently sits on disk.
#'
#' Exists so a run can state what it is actually working from. The FWC endpoints
#' behind fn.pull_seagrass_fwc() and fn.pull_reeflocations() serve the current
#' compilation, not a pinned version, and fn.pull_all() skips inputs already
#' present -- so which data you hold depends on when you first cloned. Without a
#' fingerprint two people run identical code, see identical output, and can be
#' working from different inputs with nothing to tell them apart.
#'
#' base R only, like the rest of this file, which shapes the three counts:
#'
#'   shapefile   feature count from the .shx index: 8 bytes per record after a
#'               100-byte header, so (size - 100) / 8 is exact. Verified against
#'               terra::vect() on all three shapefiles in this project.
#'   CSV         LINE count, not record count. No cheap base-R method reproduces
#'               read.csv()'s record count on files with embedded newlines --
#'               count.fields() is off by one on both artificial reef tables --
#'               and a number that is quietly wrong is worse than one that is
#'               honestly labelled. reeflocations.csv is 4,550 lines for the
#'               4,548 records FWC reports. A line count is still an exact,
#'               reproducible function of the file, which is what drift
#'               detection needs.
#'   directory   file count.
#'
#' @param row One row of fn.data_manifest().
#' @param dir.data Data root.
#' @param file.gdb Explicit geodatabase path, or NULL to discover one.
#' @return list(n, unit, bytes, md5). Any element may be NA when the input is
#'   absent or the count cannot be taken; md5 covers the primary file only.
fn.input_fingerprint <- function(row, dir.data, file.gdb = NULL) {

  none <- list(n = NA_integer_, unit = fn.input_unit(row$type),
               bytes = NA_real_, md5 = NA_character_)

  found <- fn.input_found(row, dir.data, file.gdb)
  if (!isTRUE(found$present)) return(none)

  p <- if (nzchar(row$path)) file.path(dir.data, row$path) else dir.data

  # bytes and file count over whatever the input actually is
  tree <- if (row$type == "gdb") list.files(found$found, recursive = TRUE, full.names = TRUE)
          else if (row$type == "file") found$found
          else list.files(p, recursive = TRUE, full.names = TRUE)

  bytes <- sum(file.size(tree), na.rm = TRUE)

  n <- switch(row$type,

    "dir.shp" = {
      shx <- sub("\\.shp$", ".shx", found$found)
      if (file.exists(shx)) as.integer((file.size(shx) - 100) / 8) else NA_integer_
    },

    "file" = {
      con <- file(found$found, "r")
      on.exit(close(con), add = TRUE)
      k <- 0L
      repeat {
        chunk <- readLines(con, n = 50000L, warn = FALSE)
        if (!length(chunk)) break
        k <- k + length(chunk)
      }
      k
    },

    "dir.asc" = length(list.files(p, pattern = "\\.asc$", recursive = TRUE)),

    length(tree)   # dir.any and gdb
  )

  # md5 of the primary file only - the .shp or the .csv - not the whole tree.
  # It is the thing that would actually change, and hashing 274 MB to say so
  # would not be worth the wait.
  primary <- if (row$type %in% c("dir.shp", "file")) found$found else NA_character_
  md5 <- if (!is.na(primary)) unname(tools::md5sum(primary)) else NA_character_

  list(n = as.integer(n), unit = fn.input_unit(row$type),
       bytes = as.numeric(bytes), md5 = md5)
}


#' Print the acquisition instructions for one input.
#'
#' This is what a new user sees when something is missing, so it says all four
#' things at once: what it is, where it comes from, exactly where to put it, and
#' what shape it has to be in.
#'
#' @param row One row of fn.data_manifest().
#' @param dir.data Data root.
#' @return NULL, invisibly. Called for the printed output.
fn.input_instructions <- function(row, dir.data) {

  target <- if (nzchar(row$path)) file.path(dir.data, row$path) else dir.data

  cat(sprintf("\n  %s  (section %s, %s)\n", row$key, row$section,
              if (row$required) "REQUIRED" else "optional"))
  cat(sprintf("    source : %s\n", row$url))
  cat(sprintf("    terms  : %s\n", row$terms))
  cat(sprintf("    place  : %s\n", target))
  cat(sprintf("    expects: %s\n", row$expects))
  if (!is.na(row$pull))
    cat(sprintf("    fetch  : fn.pull_all(dir.data), or %s() directly\n", row$pull))

  invisible(NULL)
}


#' Report which inputs are present and how to obtain the rest.
#'
#' Call this once at the top of the driver. A new user then learns everything
#' that is missing in one pass, rather than 40 minutes into a run.
#'
#' @param dir.data Data root.
#' @param file.gdb Explicit geodatabase path, or NULL to discover one.
#' @param stop.on.missing Stop when a REQUIRED input is missing.
#' @param fingerprint Count the items in each present input and compare against
#'   fn.reference_counts(). Set FALSE to skip the I/O and only test existence.
#' @param verbose Print the status table and the instructions.
#' @return Named logical vector, one element per manifest key, invisibly. It
#'   carries a "drift" attribute naming any input whose count has moved.
fn.check_inputs <- function(dir.data, file.gdb = NULL, stop.on.missing = TRUE,
                            fingerprint = TRUE, verbose = TRUE) {

  man <- fn.data_manifest()
  st  <- lapply(seq_len(nrow(man)), function(i)
                fn.input_found(man[i, ], dir.data, file.gdb))

  man$present <- vapply(st, function(x) isTRUE(x$present), logical(1))
  man$found   <- vapply(st, function(x) as.character(x$found), character(1))

  # What is actually on disk, against what produced the verified basemaps.
  ref <- fn.reference_counts()
  man$n.ref <- unname(ref[man$key])
  man$n     <- NA_integer_
  man$unit  <- vapply(man$type, fn.input_unit, character(1), USE.NAMES = FALSE)

  if (isTRUE(fingerprint)) {
    fp <- lapply(seq_len(nrow(man)), function(i)
                 fn.input_fingerprint(man[i, ], dir.data, file.gdb))
    man$n <- vapply(fp, function(x) as.integer(x$n), integer(1))
  }

  # DRIFT, not an error: FWC revises these layers, so a difference is often
  # legitimate. It just should not pass unnoticed.
  man$drift <- man$present & !is.na(man$n) & !is.na(man$n.ref) & man$n != man$n.ref

  if (verbose) {
    cat(sprintf("\nInput data check -- %s\n", dir.data))
    cat(strrep("-", 76), "\n", sep = "")
    cat(sprintf("  %-28s %-5s %-9s %-7s %-16s %s\n",
                "input", "sect", "need", "source", "holding", "status"))
    for (i in seq_len(nrow(man)))
      cat(sprintf("  %-28s %-5s %-9s %-7s %-16s %s\n",
                  man$key[i], man$section[i],
                  if (man$required[i]) "required" else "optional",
                  man$how[i],
                  if (is.na(man$n[i])) "-"
                    else sprintf("%s %s", format(man$n[i], big.mark = ","), man$unit[i]),
                  if (!man$present[i]) "MISSING" else if (man$drift[i]) "DRIFT" else "OK"))
    cat(strrep("-", 76), "\n", sep = "")
  }

  missing.auto <- man[!man$present & man$how == "auto", ]
  missing.man  <- man[!man$present & man$how != "auto", ]

  if (verbose && nrow(missing.auto) > 0) {
    cat(sprintf("\n%d input(s) download themselves. Run:  fn.pull_all(dir.data)\n",
                nrow(missing.auto)))
    for (i in seq_len(nrow(missing.auto)))
      fn.input_instructions(missing.auto[i, ], dir.data)
  }

  if (verbose && nrow(missing.man) > 0) {
    cat(sprintf("\n%d input(s) must be obtained by hand:\n", nrow(missing.man)))
    for (i in seq_len(nrow(missing.man)))
      fn.input_instructions(missing.man[i, ], dir.data)
  }

  short <- man$key[!man$present & man$required]
  if (length(short) > 0) {
    msg <- paste0("Missing required input(s): ", paste(short, collapse = ", "),
                  ". See the instructions above, or README.md > Getting the data.")
    if (isTRUE(stop.on.missing)) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }

  if (any(man$drift)) {
    d <- man[man$drift, ]
    warning("Input(s) differ from the reference the verified basemaps were built from:\n",
            paste(sprintf("  %-28s %s %s, expected %s",
                          d$key, format(d$n, big.mark = ","), d$unit,
                          format(d$n.ref, big.mark = ",")),
                  collapse = "\n"),
            "\nThe FWC layers are revised periodically, so this may be legitimate - ",
            "but your grids will not match a run built from the reference data. ",
            "See README.md > Getting the data, and data/PROVENANCE.tsv for what was ",
            "fetched when.", call. = FALSE)
  }

  if (verbose) {
    opt <- man$key[!man$present & !man$required]
    if (length(opt) > 0)
      cat(sprintf("\nProceeding without optional input(s): %s.\nThe affected sections degrade rather than fail -- see README.md > Getting the data.\n",
                  paste(opt, collapse = ", ")))
    else if (all(man$present) && !any(man$drift))
      cat("\nAll inputs present, and matching the reference counts.\n")
  }

  out <- stats::setNames(man$present, man$key)
  attr(out, "drift") <- man$key[man$drift]
  invisible(out)
}


#' Append one downloaded input to the provenance record.
#'
#' data/PROVENANCE.tsv accumulates a row per successful download -- when, from
#' where, and a fingerprint of what arrived. It is append-only, so it is the
#' history of this working copy rather than a snapshot.
#'
#' It lives under data/ and is therefore gitignored, deliberately: it describes
#' YOUR copy, so committing it would mean merge conflicts over machine-specific
#' facts. The shared reference is fn.reference_counts(), which is tracked.
#'
#' @param row One row of fn.data_manifest().
#' @param dir.data Data root.
#' @return Path to the provenance file, invisibly.
fn.record_provenance <- function(row, dir.data) {

  f  <- file.path(dir.data, "PROVENANCE.tsv")
  fp <- fn.input_fingerprint(row, dir.data)

  if (!file.exists(f))
    cat("fetched\tinput\tn\tunit\tbytes\tmd5\turl\n", file = f)

  cat(sprintf("%s\t%s\t%s\t%s\t%.0f\t%s\t%s\n",
              format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), row$key,
              if (is.na(fp$n)) "NA" else fp$n, fp$unit,
              if (is.na(fp$bytes)) 0 else fp$bytes,
              if (is.na(fp$md5)) "NA" else fp$md5, row$url),
      file = f, append = TRUE)

  invisible(f)
}


#' Run every automatic download that is still needed.
#'
#' Skips anything already on disk, so it is safe to re-run. Each download is
#' wrapped so one unreachable endpoint does not abort the rest -- a failure is
#' reported with the page to fetch that dataset from by hand. Whatever is
#' fetched is appended to data/PROVENANCE.tsv.
#'
#' @param dir.data Data root.
#' @param overwrite Re-download inputs that are already present.
#' @return Named logical vector of per-dataset success, invisibly.
fn.pull_all <- function(dir.data, overwrite = FALSE) {

  man <- fn.data_manifest()
  man <- man[man$how == "auto", ]
  ok  <- stats::setNames(logical(nrow(man)), man$key)

  for (i in seq_len(nrow(man))) {

    row <- man[i, ]

    if (!overwrite && isTRUE(fn.input_found(row, dir.data)$present)) {
      message(sprintf("%-24s already present, skipping", row$key))
      ok[i] <- TRUE
      next
    }

    dir.target <- file.path(dir.data, row$pull.dir)
    if (!dir.exists(dir.target)) dir.create(dir.target, recursive = TRUE)

    message(sprintf("%-24s downloading...", row$key))
    res <- try(do.call(row$pull, list(dir.target)), silent = TRUE)

    if (inherits(res, "try-error")) {
      warning(sprintf("%s failed to download: %s\n  Obtain it by hand from %s",
                      row$key, conditionMessage(attr(res, "condition")), row$url),
              call. = FALSE)
    } else {
      ok[i] <- isTRUE(fn.input_found(row, dir.data)$present)
      if (!ok[i]) {
        warning(sprintf("%s downloaded but does not match the expected shape: %s",
                        row$key, row$expects), call. = FALSE)
      } else {
        fn.record_provenance(row, dir.data)
      }
    }
  }

  invisible(ok)
}
