import calliope
model = calliope.Model("model.yaml", scenario="only_gbr") # or however you remove all non-UK nodes
techs = model.inputs.name.to_series()
print(techs)