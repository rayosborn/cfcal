from __future__ import absolute_import
from . import define_model

def plugin_menu(parent):
    menu = 'CF'
    actions = [('Set CF Parameters', define_model.show_dialog))
    return menu, actions
