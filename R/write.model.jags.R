write.model.jags <-
  function (model, name_model, intitled, Data, param = "value") {
    model <- replace.inprod(body(model), name_model, Data, param)
    writeLines(model, intitled)
  }
