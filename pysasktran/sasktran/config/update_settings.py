import sys
from . import load_user_config, save_user_config


def update_glossac_directory(fullfilename: str):
    user_config = load_user_config()
    user_config['glossac_file'] = fullfilename
    save_user_config(user_config)
    print('Updated Glosacc file to <{:s}>'.format(fullfilename))


if __name__ == "__main__":
    args = sys.argv
    if (len(args) > 1):
        command = args[1].lower()
        if command == 'set_glossac_file':
            filename = args[2]
            update_glossac_directory(filename)
        else:
            print('Unrecognized command')
            exit(1)
