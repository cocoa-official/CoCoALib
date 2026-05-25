import json
import os

def main():
    directories = ['packages', 'tests', 'CoCoAManual']
    paths = []

    for d in directories:
        if os.path.exists(d):
            for root, _, files in os.walk(d):
                for file in files:
                    app_path = os.path.join(root, file).replace('\\', '/')
                    paths.append(app_path)

    try:
        with open('cocoa-fs.json', 'w', encoding='utf-8') as f:
            json.dump(paths, f, separators=(',', ':'))
        print(f"Successfully wrote {len(paths)} files to cocoa-fs.json")
    except Exception as e:
        print(f"Failed to write cocoa-fs.json: {e}")

if __name__ == "__main__":
    main()