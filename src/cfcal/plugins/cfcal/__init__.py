from . import define_model, fit_model

def plugin_menu():
    menu = 'CF'
    actions = [('Set CF Parameters', define_model.show_dialog),
               ('Fit CF Model', fit_model.show_dialog)]
    return menu, actions
