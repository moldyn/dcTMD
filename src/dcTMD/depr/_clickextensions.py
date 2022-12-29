# -*- coding: utf-8 -*-
"""
Define new click class.

Mutex by Stephen Rauch, omni: https://stackoverflow.com/questions/44247099
"""


import click


class Mutex(click.Option):
    """Add mutually exclusive click options."""

    def __init__(self, *args, **kwargs):
        self.not_required_if: list = kwargs.pop('not_required_if')

        assert self.not_required_if, "'not_required_if' parameter required"
        kwargs['help'] = (
            kwargs.get('help', '') + ', option is mutually exclusive with ' +
            ', '.join(self.not_required_if)
        ).strip()
        super(Mutex, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        current_opt: bool = self.name in opts
        for mutex_opt in self.not_required_if:
            if mutex_opt in opts:
                if current_opt:
                    raise click.UsageError(
                        "Illegal usage: '" + str(self.name) +
                        "' is mutually exclusive with '" + str(mutex_opt) +
                        "'."
                    )
                else:
                    self.prompt = None
        return super(Mutex, self).handle_parse_result(ctx, opts, args)
