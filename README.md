# Sylens' Blog

Source code for my blog, use hexo to generate static pages.

The theme is [Volantis](https://volantis.js.org/), a very powerful theme with many additional function.

## Structure

To fullfill the need for multi-language support, which is difficut for hexo, I referenced [the method from a blog](https://loli.fj.cn/2023/09/09/Hexo%E9%85%8D%E7%BD%AE%E5%A4%9A%E8%AF%AD%E8%A8%80/#%E9%85%8D%E7%BD%AE%E7%BF%BB%E8%AF%91%E7%9A%84%E9%A1%B5%E9%9D%A2%EF%BC%88%E5%8F%AF%E9%80%89%EF%BC%89), that is, to build 2 seperate blog project, and set different languange for them, then merge the compile content together:

```bash
.
├── blog_cn # Chinese blog
├── blog_en  # English blog
├── pixi.lock
├── pixi.toml # pixi dep and task for deploy
├── aider.sh # script use aider to translate zh posts to en ones
├── .devcontainer # config files for using DevPod to build workspace
└── README.md
```

[pixi](http://pixi.sh/) was used to deploy needed nodejs dependency and run build/publish tasks

### Patch

Add a patch file to fix some small issue of `hexo-theme-volantis`:
- Invalid `app.js` path when `:root` has value
- Unused code in the visitor statistics section at the bottom of the page
- Incomplete English localization caused by Hexo framework itself.

## Deploy & Develop

### Using Pixi

- Install [pixi](http://pixi.sh/)
- fork this repo, enther the directory, and type `pixi run build`

### DevPod

This repo no has DevContainer config files, install [DevPod](https://devpod.sh/) Desktop or Cli, and then specify the path of this repo can start workspace for this repo